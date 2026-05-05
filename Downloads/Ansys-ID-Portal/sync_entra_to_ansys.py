"""
Entra ID → Ansys ID Portal Sync
---------------------------------
Synchronizes members of an Azure Entra ID (AAD) security group to either:
  - An Ansys ID Portal Account (AccountUser membership), or
  - An Ansys ID Portal Group (GroupMember membership).

Usage:
    python sync_entra_to_ansys.py ^
        --entra-domain "contoso.onmicrosoft.com" ^
        --entra-group "MySecurityGroup" ^
        --entra-client-id "your-app-client-id" ^
        --entra-client-secret "your-app-secret" ^
        [--ansys-pat "eyJ..."] ^
        --target-type account ^
        --account-number "my_account_number" ^
        [--group-name "MyAnsysGroup"] ^
        [--dry-run] ^
        [--send-email] ^
        [--log-file "C:\\Logs\\ansys_sync.log"] ^
        [--smtp-relay "smtp.contoso.com"] ^
        [--smtp-port 25] ^
        [--smtp-from "ansys-sync@contoso.com"] ^
        [--log-recipient "it-admin@contoso.com"] ^
        [--warn-recipient "it-alerts@contoso.com"]

    (Use backslash line continuation on Linux/macOS instead of ^)

    --entra-client-id and --entra-client-secret may be omitted from the
    command line and the script will prompt for them interactively.
    The secret prompt does not echo to the screen.

Arguments:
    --entra-domain        Entra ID tenant domain, e.g. contoso.onmicrosoft.com
    --entra-group         Display name of the Entra ID security group
    --entra-client-id     App registration client ID for Graph API access
                          (prompted if not supplied)
    --entra-client-secret App registration client secret for Graph API access
                          (prompted silently if not supplied)
    --ansys-pat           Ansys ID Portal PAT for B2C authentication
                          (prompted interactively if not supplied)
    --target-type         Either "account" or "group"
    --account-number      Ansys ID Portal account number (must already exist)
    --group-name          Ansys group name (required when --target-type is group)
    --dry-run             Print planned changes without applying them
    --send-email          Send Ansys ID invitation emails when adding users
                          (emails are suppressed by default)
    --reset-credentials   Clear all cached credentials and prompt for new ones

    Logging / Notification:
    --log-file            Path to append the full run log to
    --smtp-relay          Hostname or IP of the SMTP relay/gateway
    --smtp-port           SMTP port (default: 25)
    --smtp-from           Envelope From address for all outbound emails
    --log-recipient       Comma-separated address(es) for the full log email
    --warn-recipient      Comma-separated address(es) for warnings/errors email

No environment variables are required. All credentials are either passed
as CLI arguments or prompted interactively at runtime.

Dependencies:
    pip install msal httpx requests

The Entra ID app registration needs the following Graph API application
permissions:
    Group.Read.All
    User.Read.All

API Reference: https://iam.ansys.com/swagger/AnsysId/swagger.json
"""

import sys
import json
import asyncio
import argparse
import getpass
import logging
import smtplib
from datetime import datetime, timezone
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from io import StringIO
from pathlib import Path
from urllib.parse import urlparse, parse_qs

import os

import requests
import httpx
import msal

try:
    from cryptography.hazmat.primitives.ciphers.aead import AESGCM
    from cryptography.hazmat.primitives.kdf.pbkdf2 import PBKDF2HMAC
    from cryptography.hazmat.primitives import hashes
    CRYPTO_AVAILABLE = True
except ImportError:
    CRYPTO_AVAILABLE = False

# ============================================================
# GLOBALS — populated after arg parsing
# ============================================================
log = logging.getLogger("ansys_sync")

# Captures WARNING+ records for the warnings-only email.
_warn_records: list[logging.LogRecord] = []

# Structured sync summary used in email subjects and bodies.
_sync_summary: dict = {
    "started_at":  None,
    "finished_at": None,
    "entra_group": None,
    "target_type": None,
    "target_name": None,
    "added":       [],
    "removed":     [],
    "skipped":     [],   # members with no primary email (mail) attribute
    "errors":      [],
    "dry_run":     False,
    "exit_status": "success",
}

# ============================================================
# LOGGING SETUP
# ============================================================

class _WarnCollector(logging.Handler):
    """Captures WARNING+ records into _warn_records for the alert email."""
    def emit(self, record: logging.LogRecord) -> None:
        if record.levelno >= logging.WARNING:
            _warn_records.append(record)


def setup_logging(log_file: str | None) -> StringIO:
    """
    Attaches four handlers to the logger:
      - Console       INFO+   (always)
      - StringIO      DEBUG+  (always — captured for log email)
      - File          DEBUG+  (only if --log-file supplied)
      - WarnCollector WARNING+ (always — captured for warn email)

    Returns the StringIO buffer so main() can attach it to emails.
    """
    fmt = logging.Formatter(
        "%(asctime)s [%(levelname)-8s] %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )
    log.setLevel(logging.DEBUG)

    # Console
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)
    log.addHandler(ch)

    # In-memory buffer
    buf = StringIO()
    sh  = logging.StreamHandler(buf)
    sh.setLevel(logging.DEBUG)
    sh.setFormatter(fmt)
    log.addHandler(sh)

    # File (appended across runs)
    if log_file:
        fh = logging.FileHandler(log_file, encoding="utf-8")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(fmt)
        log.addHandler(fh)
        log.info("Log file: %s", log_file)

    # Warning collector
    log.addHandler(_WarnCollector())

    return buf


# ============================================================
# SMTP HELPERS
# ============================================================

def _send_email(smtp_relay: str, smtp_port: int, from_addr: str,
                to_addrs: list[str], subject: str, body: str) -> None:
    """Send a plain-text email via an unauthenticated SMTP relay."""
    if not smtp_relay or not to_addrs or not from_addr:
        return
    msg            = MIMEMultipart("alternative")
    msg["Subject"] = subject
    msg["From"]    = from_addr
    msg["To"]      = ", ".join(to_addrs)
    msg.attach(MIMEText(body, "plain", "utf-8"))
    try:
        with smtplib.SMTP(smtp_relay, smtp_port, timeout=15) as server:
            server.sendmail(from_addr, to_addrs, msg.as_string())
        log.info("Email sent to %s — subject: %s", to_addrs, subject)
    except Exception as exc:
        # Never let an email failure mask the sync result
        log.warning("Failed to send email to %s: %s", to_addrs, exc)


def send_log_email(smtp_relay: str, smtp_port: int, from_addr: str,
                   recipients: list[str], log_buffer: StringIO) -> None:
    """Send the full run log to --log-recipient addresses."""
    if not recipients:
        return
    s   = _sync_summary
    tag = "[DRY-RUN] " if s["dry_run"] else ""
    subject = (
        f"{tag}Ansys Sync {s['exit_status'].upper()} — "
        f"{s['entra_group']} → {s['target_type']}:{s['target_name']} — "
        f"{datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}"
    )
    body = "\n".join([
        "Ansys ID Portal Sync — Full Log",
        "=" * 60,
        f"Started :  {s['started_at']}",
        f"Finished:  {s['finished_at']}",
        f"Status  :  {s['exit_status'].upper()}"
            + (" (DRY-RUN)" if s["dry_run"] else ""),
        f"Source  :  Entra group '{s['entra_group']}'",
        f"Target  :  {s['target_type']} '{s['target_name']}'",
        "",
        f"Added   : {len(s['added'])}  {sorted(s['added'])}",
        f"Removed : {len(s['removed'])}  {sorted(s['removed'])}",
        f"Skipped : {len(s['skipped'])}  {sorted(s['skipped'])}",
        f"Errors  : {len(s['errors'])}",
        "",
        "=" * 60,
        "FULL LOG",
        "=" * 60,
        log_buffer.getvalue(),
    ])
    _send_email(smtp_relay, smtp_port, from_addr, recipients, subject, body)


def send_warn_email(smtp_relay: str, smtp_port: int, from_addr: str,
                    recipients: list[str]) -> None:
    """
    Send a warnings/errors summary to --warn-recipient addresses.
    Only sent when there are actual warnings or errors to report.
    """
    if not recipients:
        return
    if not _warn_records and not _sync_summary["errors"]:
        log.info("No warnings or errors — skipping alert email.")
        return

    s   = _sync_summary
    tag = "[DRY-RUN] " if s["dry_run"] else ""
    subject = (
        f"{tag}⚠ Ansys Sync WARNINGS — "
        f"{s['entra_group']} → {s['target_type']}:{s['target_name']} — "
        f"{datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}"
    )
    lines = [
        "Ansys ID Portal Sync — Warnings & Errors Summary",
        "=" * 60,
        f"Run time : {s['started_at']} → {s['finished_at']}",
        f"Status   : {s['exit_status'].upper()}"
            + (" (DRY-RUN)" if s["dry_run"] else ""),
        f"Source   : Entra group '{s['entra_group']}'",
        f"Target   : {s['target_type']} '{s['target_name']}'",
        "",
    ]
    if s["skipped"]:
        lines.append(
            f"SKIPPED USERS ({len(s['skipped'])}) "
            "— no primary email (mail) attribute in Entra:"
        )
        for entry in sorted(s["skipped"]):
            lines.append(f"  • {entry}")
        lines.append("")
    if s["errors"]:
        lines.append(f"ERRORS ({len(s['errors'])}):")
        for err in s["errors"]:
            lines.append(f"  • {err}")
        lines.append("")
    if _warn_records:
        fmt = logging.Formatter(
            "%(asctime)s [%(levelname)s] %(message)s",
            datefmt="%Y-%m-%dT%H:%M:%S",
        )
        lines.append(f"LOG WARNINGS/ERRORS ({len(_warn_records)}):")
        for r in _warn_records:
            lines.append(f"  {fmt.format(r)}")

    _send_email(smtp_relay, smtp_port, from_addr, recipients,
                subject, "\n".join(lines))


# ============================================================
# ANSYS B2C CONFIGURATION
# ============================================================
POLICY       = "B2C_1A_ANSYSID_SIGNUP_SIGNIN"
REDIRECT_URI = "https://login.microsoftonline.com/common/oauth2/nativeclient"

# Production endpoints and credentials -- all hardcoded.
AADB2C_ENV       = "ansysaccount"
TENANT_NAME      = "ansysaccount.onmicrosoft.com"
ANSYS_IAM_BASE   = "https://iam.ansys.com"
TOKEN_CACHE_FILE = Path(".token_cache.json")
CLIENT_ID        = "7ef9d43e-407a-40a2-94d7-14eb40416af8"
ANSYS_SCOPE      = "https://ansysaccount.onmicrosoft.com/AnsysID/Authentication"
GRAPH_BASE       = "https://graph.microsoft.com/v1.0"

# ============================================================
# ENTRA ID CREDENTIAL RESOLUTION
# ============================================================

# ============================================================
# CREDENTIAL CACHE — AES-256-GCM encrypted file, cross-platform
# ============================================================
# Works identically on Windows, macOS, and headless Linux.
# Credentials are stored in an encrypted JSON file protected by a master
# password that is prompted once per run and never written to disk.
#
# File layout (all binary fields base64-encoded, stored as JSON):
#   { "salt": <16 bytes>, "nonce": <12 bytes>, "ciphertext": <bytes> }
#
# Key derivation: PBKDF2-SHA256, 600,000 iterations, 32-byte output.
# Encryption:     AES-256-GCM (authenticated — detects tampering).
#
# Dependency: pip install cryptography
# ============================================================

CRED_CACHE_FILE = Path(".credential_cache.json")

# Master password held in memory for the duration of the run only.
_master_password: str | None = None


def _get_master_password(confirm: bool = False) -> str:
    """
    Return the master cache password, prompting if not yet set this run.
    If confirm=True, prompts twice (used when creating a new cache file).
    """
    global _master_password
    if _master_password:
        return _master_password
    while True:
        pw = getpass.getpass(
            "Enter credential cache password"
            " (protects cached Entra and Ansys credentials): "
        )
        if not pw:
            log.warning("Password cannot be empty.")
            continue
        if confirm:
            pw2 = getpass.getpass("Confirm credential cache password: ")
            if pw != pw2:
                log.warning("Passwords do not match — try again.")
                continue
        _master_password = pw
        return pw


def _derive_key(password: str, salt: bytes) -> bytes:
    """Derive a 32-byte AES key from password + salt using PBKDF2-SHA256."""
    kdf = PBKDF2HMAC(
        algorithm=hashes.SHA256(),
        length=32,
        salt=salt,
        iterations=600_000,
    )
    return kdf.derive(password.encode("utf-8"))


def _load_cred_cache() -> dict:
    """
    Load and decrypt the credential cache file.
    Returns {} if the file doesn't exist, the library is unavailable,
    or decryption fails (wrong password or tampered file).
    """
    if not CRYPTO_AVAILABLE:
        log.debug("cryptography not installed — credential caching disabled.")
        return {}
    if not CRED_CACHE_FILE.exists():
        return {}
    try:
        import base64
        raw   = json.loads(CRED_CACHE_FILE.read_text())
        salt  = base64.b64decode(raw["salt"])
        nonce = base64.b64decode(raw["nonce"])
        ct    = base64.b64decode(raw["ciphertext"])
        key   = _derive_key(_get_master_password(), salt)
        pt    = AESGCM(key).decrypt(nonce, ct, None)
        return json.loads(pt.decode("utf-8"))
    except Exception as exc:
        log.warning(
            "Could not decrypt credential cache (%s). "
            "Credentials will be re-prompted.", exc,
        )
        global _master_password
        _master_password = None   # allow re-entry of password
        return {}


def _save_cred_cache(creds: dict) -> None:
    """Encrypt and write the credential cache file."""
    if not CRYPTO_AVAILABLE:
        return
    try:
        import base64
        is_new = not CRED_CACHE_FILE.exists()
        salt   = os.urandom(16)
        nonce  = os.urandom(12)
        key    = _derive_key(_get_master_password(confirm=is_new), salt)
        ct     = AESGCM(key).encrypt(nonce, json.dumps(creds).encode(), None)
        CRED_CACHE_FILE.write_text(json.dumps({
            "salt":       base64.b64encode(salt).decode(),
            "nonce":      base64.b64encode(nonce).decode(),
            "ciphertext": base64.b64encode(ct).decode(),
        }))
        log.debug("Credential cache written to %s", CRED_CACHE_FILE)
    except Exception as exc:
        log.warning("Could not write credential cache: %s", exc)


def _cred_get(key: str) -> str | None:
    """Retrieve a single credential from the encrypted cache."""
    return _load_cred_cache().get(key)


def _cred_set(key: str, value: str) -> None:
    """Store a single credential in the encrypted cache."""
    creds = _load_cred_cache()
    creds[key] = value
    _save_cred_cache(creds)


def reset_cached_credentials() -> None:
    """Delete the credential cache file and clear the in-memory password."""
    if CRED_CACHE_FILE.exists():
        CRED_CACHE_FILE.unlink()
        log.info("Credential cache file deleted: %s", CRED_CACHE_FILE)
    else:
        log.info("No credential cache file found.")
    global _master_password
    _master_password = None


def resolve_entra_credentials(args: argparse.Namespace) -> tuple[str, str]:
    """
    Return (client_id, client_secret).

    Resolution order for each value:
      1. CLI argument (--entra-client-id / --entra-client-secret)
      2. OS keychain cache (Credential Manager on Windows)
      3. Interactive prompt (stored to keychain for future runs)
    """
    # ── Client ID ───────────────────────────────────────────────────────────
    client_id = args.entra_client_id or _cred_get("entra_client_id")
    if not client_id:
        client_id = input(
            "Enter Entra App Registration Client ID: "
        ).strip()
        if not client_id:
            raise RuntimeError("Entra client ID is required.")
        _cred_set("entra_client_id", client_id)
        log.info("Entra client ID saved to credential cache.")
    else:
        log.debug("Entra client ID retrieved.")

    # ── Client Secret ────────────────────────────────────────────────────────
    client_secret = args.entra_client_secret or _cred_get("entra_client_secret")
    if not client_secret:
        client_secret = getpass.getpass(
            "Enter Entra App Registration Client Secret: "
        ).strip()
        if not client_secret:
            raise RuntimeError("Entra client secret is required.")
        _cred_set("entra_client_secret", client_secret)
        log.info("Entra client secret saved to credential cache.")
    else:
        log.debug("Entra client secret retrieved.")

    return client_id, client_secret.strip()


# ============================================================
# ENVIRONMENT CONFIGURATION
# ============================================================

def configure_environment() -> None:
    """Log the fixed production endpoint configuration."""
    log.info(
        "B2C: %s  ID Portal: %s  Scope: %s",
        f"https://{AADB2C_ENV}.b2clogin.com", ANSYS_IAM_BASE, ANSYS_SCOPE,
    )


# ============================================================
# TOKEN CACHE HELPERS
# ============================================================

def load_token_cache() -> msal.SerializableTokenCache:
    cache = msal.SerializableTokenCache()
    if TOKEN_CACHE_FILE.exists():
        cache.deserialize(TOKEN_CACHE_FILE.read_text())
    return cache


def save_token_cache(cache: msal.SerializableTokenCache) -> None:
    if cache.has_state_changed:
        TOKEN_CACHE_FILE.write_text(cache.serialize())


def extract_email_from_token(token: str) -> str | None:
    """
    Extract the email claim from a JWT id_token or access_token without
    requiring any external library. JWTs are three base64url-encoded
    segments separated by dots — the middle segment is the payload.
    """
    import base64
    try:
        payload_segment = token.split(".")[1]
        # base64url requires padding to a multiple of 4
        padding = 4 - len(payload_segment) % 4
        payload_segment += "=" * (padding % 4)
        payload = json.loads(
            base64.urlsafe_b64decode(payload_segment).decode("utf-8")
        )
        email = payload.get("email") or payload.get("preferred_username")
        if email:
            log.debug("Extracted email from token: %s", email)
        return email
    except Exception as exc:
        log.debug("Could not extract email from token: %s", exc)
        return None


# ============================================================
# ANSYS B2C AUTH
# ============================================================

def build_authorization_url(pat: str, client_id: str,
                            resource_scope: str) -> str:
    scope = f"openid offline_access {resource_scope}"
    return (
        f"https://{AADB2C_ENV}.b2clogin.com/"
        f"{TENANT_NAME}/oauth2/v2.0/authorize"
        f"?p={POLICY}"
        f"&client_id={client_id}"
        f"&nonce=defaultNonce"
        f"&redirect_uri={REDIRECT_URI}"
        f"&scope={scope}"
        f"&response_type=code"
        f"&id_token_hint={pat}"
    )


async def _get_auth_code(client: httpx.AsyncClient, auth_url: str) -> str:
    resp     = await client.get(auth_url, follow_redirects=False)
    location = resp.headers.get("Location")
    if not location:
        # B2C returned an HTML error page instead of redirecting.
        # Extract the Detail field from the embedded GLOBALEX JSON object
        # which contains the specific AADB2C error code and description.
        detail = "unknown — could not parse B2C error detail"
        try:
            import re as _re
            match = _re.search(
                r'"Detail"\s*:\s*"([^"]+)"',
                resp.text,
            )
            if match:
                detail = match.group(1)
        except Exception:
            pass
        log.error(
            "B2C authorization endpoint returned HTTP %s instead of "
            "redirecting. B2C detail: %s",
            resp.status_code,
            detail,
        )
        raise RuntimeError(
            f"B2C rejected the authorization request: {detail}"
        )
    parsed = urlparse(location)
    params = parse_qs(parsed.query)
    if "code" in params:
        return params["code"][0]
    if "error" in params:
        error       = params.get("error", ["unknown"])[0]
        description = params.get("error_description", ["No description"])[0]
        raise RuntimeError(f"B2C error '{error}': {description}")
    raise RuntimeError(
        f"Redirect received but no code or error found. Location: {location}"
    )


async def get_ansys_access_token(client_id: str,
                                resource_scope: str,
                                pat_arg: str | None = None) -> tuple[str, str | None]:
    """
    Return (token, caller_email).

    token        — access_token if issued, otherwise id_token as fallback.
    caller_email — email extracted from the token JWT payload, used to
                   call GET /User/Details for account UUID resolution.

    MSAL manages openid and offline_access internally — they must NOT be
    included in the scopes list passed to MSAL methods, but they ARE
    included explicitly in the authorization URL scope parameter so that
    B2C knows to issue an ID token and refresh token.
    """
    # Scope list for MSAL calls — resource scope only.
    # MSAL implicitly adds openid and offline_access.
    msal_scopes = [resource_scope]

    cache     = load_token_cache()
    authority = f"https://{AADB2C_ENV}.b2clogin.com/{TENANT_NAME}/{POLICY}"
    app       = msal.PublicClientApplication(
        client_id=client_id, authority=authority, token_cache=cache
    )

    accounts = app.get_accounts()
    if accounts:
        result = app.acquire_token_silent(msal_scopes, account=accounts[0])
        if result and "access_token" in result:
            save_token_cache(cache)
            log.info("Ansys token retrieved from cache.")
            email = (extract_email_from_token(result.get("id_token", ""))
                     or extract_email_from_token(result["access_token"]))
            return result["access_token"], email
        if result and "id_token" in result and not result.get("error"):
            save_token_cache(cache)
            log.info("Ansys id_token retrieved from cache.")
            email = extract_email_from_token(result["id_token"])
            return result["id_token"], email

    log.warning("No valid Ansys token cached. PAT authentication required.")
    # Resolution order: CLI arg → credential cache → interactive prompt
    if pat_arg:
        pat = pat_arg
        log.info("Ansys PAT supplied via command line argument.")
    else:
        pat = _cred_get("ansys_pat")
        if pat:
            log.info("Ansys PAT retrieved from credential cache.")
        else:
            log.warning(
                "The Ansys PAT is a JWT token obtained from the Ansys ID "
                "Portal — it is NOT your Ansys account password."
            )
            pat = getpass.getpass("Enter your Ansys PAT: ")
            if not pat:
                raise RuntimeError("Ansys PAT is required.")
            _cred_set("ansys_pat", pat)
            log.info("Ansys PAT saved to credential cache.")

    # build_authorization_url includes openid and offline_access explicitly
    # in the URL scope so B2C issues the correct token types.
    auth_url = build_authorization_url(pat, client_id, resource_scope)
    async with httpx.AsyncClient(timeout=15) as client:
        try:
            code = await _get_auth_code(client, auth_url)
            log.info("Authorization code retrieved.")
        except Exception as exc:
            err = f"Failed to obtain authorization code: {exc}"
            log.error(err)
            _sync_summary["errors"].append(err)
            _sync_summary["exit_status"] = "failure"
            sys.exit(1)

    result = app.acquire_token_by_authorization_code(
        code=code, scopes=msal_scopes, redirect_uri=REDIRECT_URI
    )

    if "access_token" in result:
        save_token_cache(cache)
        log.info("Ansys access token acquired and cached.")
        email = (extract_email_from_token(result.get("id_token", ""))
                 or extract_email_from_token(result["access_token"]))
        return result["access_token"], email

    # B2C may return only an id_token when the resource scope is not
    # registered in the B2C tenant (scope comes back empty). In this case
    # the id_token itself is used as the bearer credential for the API.
    if "id_token" in result and not result.get("error"):
        save_token_cache(cache)
        log.info(
            "No access token in response (scope: '%s') — "
            "using id_token as bearer credential.",
            result.get("scope", ""),
        )
        return result["id_token"]

    # Genuine failure — log every field for diagnosis
    log.error("Ansys token exchange failed. Full MSAL response:")
    for key, value in result.items():
        log.error("  %s: %s", key, value)

    err = (
        f"Ansys authentication failed — "
        f"error: {result.get('error')}  "
        f"suberror: {result.get('error_subcode') or result.get('suberror')}  "
        f"description: {result.get('error_description')}  "
        f"correlation_id: {result.get('correlation_id')}"
    )
    log.error(err)
    _sync_summary["errors"].append(err)
    _sync_summary["exit_status"] = "failure"
    sys.exit(1)


# ============================================================
# ENTRA ID — GRAPH API
# ============================================================

def _graph_headers(entra_domain: str, client_id: str,
                   client_secret: str) -> dict:
    """Obtain a Graph API token via client-credentials flow."""
    token_url = (
        f"https://login.microsoftonline.com/{entra_domain}"
        f"/oauth2/v2.0/token"
    )
    resp = requests.post(token_url, data={
        "grant_type":    "client_credentials",
        "client_id":     client_id,
        "client_secret": client_secret,
        "scope":         "https://graph.microsoft.com/.default",
    })
    resp.raise_for_status()
    return {"Authorization": f"Bearer {resp.json()['access_token']}"}


def entra_email(member: dict) -> str | None:
    """
    Returns the primary SMTP address (mail attribute) only.
    UPN is NOT used as a fallback. Members without a mail address
    are skipped with a warning and recorded in _sync_summary["skipped"].
    """
    email = member.get("mail")
    if not email:
        display = member.get("displayName", "<unknown>")
        uid     = member.get("id", "<unknown>")
        log.warning(
            "Skipping member '%s' (id: %s) — no primary email (mail) "
            "attribute set in Entra. Ensure all group members have a "
            "mail address populated.",
            display, uid,
        )
        _sync_summary["skipped"].append(f"{display} ({uid})")
        return None
    return email.lower()


def get_entra_group_members(entra_domain: str, group_name: str,
                             client_id: str,
                             client_secret: str) -> list[dict]:
    """
    Return members of the named Entra ID security group that have a
    primary email address. Uses transitiveMembers to flatten nested groups.
    """
    log.info(
        "Fetching Entra ID group '%s' in domain '%s'",
        group_name, entra_domain,
    )
    headers = _graph_headers(entra_domain, client_id, client_secret)

    # Locate the group by display name
    filter_url = (
        f"{GRAPH_BASE}/groups"
        f"?$filter=displayName eq '{group_name}' and securityEnabled eq true"
        f"&$select=id,displayName"
    )
    resp = requests.get(filter_url, headers=headers)
    resp.raise_for_status()
    groups = resp.json().get("value", [])
    if not groups:
        raise RuntimeError(
            f"Entra ID security group '{group_name}' not found."
        )
    if len(groups) > 1:
        log.warning(
            "Multiple groups matched '%s'; using first result.", group_name
        )
    group_id = groups[0]["id"]
    log.info("Found Entra group id: %s", group_id)

    # Page through all transitive (nested) members
    members: list[dict] = []
    next_url = (
        f"{GRAPH_BASE}/groups/{group_id}/transitiveMembers"
        "?$select=displayName,mail,userPrincipalName,id"
        "&$top=999"
    )
    while next_url:
        resp     = requests.get(next_url, headers=headers)
        resp.raise_for_status()
        data     = resp.json()
        members.extend(
            m for m in data.get("value", [])
            if m.get("@odata.type") == "#microsoft.graph.user"
        )
        next_url = data.get("@odata.nextLink")

    # Filter to only members with a primary email; entra_email() logs skips
    valid   = [m for m in members if entra_email(m)]
    skipped = len(members) - len(valid)

    if skipped:
        log.warning(
            "%d member(s) skipped due to missing mail attribute. "
            "They will not be synced.", skipped,
        )
    if not valid:
        raise RuntimeError(
            f"No members with a primary email address found in Entra group "
            f"'{group_name}'. Sync aborted."
        )

    log.info(
        "Retrieved %d eligible member(s) from Entra group (%d skipped).",
        len(valid), skipped,
    )
    return valid


# ============================================================
# ANSYS ID PORTAL CLIENT
# ============================================================

class AnsysIAMClient:
    """Wrapper around the Ansys ID Portal REST API."""

    def __init__(self, access_token: str, dry_run: bool = False,
                 send_email: bool = False,
                 preserve_emails: set | None = None):
        self.headers = {
            "Authorization": f"Bearer {access_token}",
            "Content-Type":  "application/json",
            "Accept":        "application/json",
        }
        self.dry_run        = dry_run
        self.send_email = send_email
        self.preserve_emails = preserve_emails or set()
        if self.preserve_emails:
            log.info("Preserved emails (never removed): %s",
                     sorted(self.preserve_emails))

    def _get(self, path: str, params: dict = None):
        url  = f"{ANSYS_IAM_BASE}{path}"
        resp = requests.get(url, headers=self.headers, params=params)
        resp.raise_for_status()
        return resp.json() if resp.content else None

    def _post(self, path: str, body=None, params: dict = None):
        url  = f"{ANSYS_IAM_BASE}{path}"
        # Send no body at all when body is None — some endpoints return 500
        # if they receive an unexpected empty JSON object.
        if body is None:
            resp = requests.post(url, headers=self.headers, params=params)
        else:
            resp = requests.post(url, headers=self.headers,
                                 json=body, params=params)
        resp.raise_for_status()
        return resp.json() if resp.content else None

    def _put(self, path: str, body=None, params: dict = None):
        url = f"{ANSYS_IAM_BASE}{path}"
        if self.dry_run:
            log.info(
                "[DRY-RUN] PUT %s  params=%s  body=%s",
                url, params,
                json.dumps(body) if body is not None else None,
            )
            return {}
        resp = requests.put(url, headers=self.headers,
                            json=body, params=params)
        resp.raise_for_status()
        return resp.json() if resp.content else None

    def _delete(self, path: str, body=None, params: dict = None):
        url = f"{ANSYS_IAM_BASE}{path}"
        if self.dry_run:
            log.info(
                "[DRY-RUN] DELETE %s  body=%s", url,
                json.dumps(body) if body is not None else None,
            )
            return None
        resp = requests.delete(url, headers=self.headers,
                               json=body, params=params)
        resp.raise_for_status()
        return resp.json() if resp.content else None

    # ── Account lookup ──────────────────────────────────────────────────────

    def find_account_uuid(self, account_number: str,
                          caller_email: str | None = None) -> str:
        """
        Resolve an account number to its UUID via GET /User/Details?email=...

        Uses the authenticated user's own account list — accessible to any
        account member without requiring an admin role. Requires the caller
        email extracted from the token JWT to identify the user.

        Aborts if the account cannot be found, which indicates the PAT
        owner is not a member of the specified account.
        """
        if not caller_email:
            err = (
                "Cannot resolve account UUID — caller email could not be "
                "extracted from the token. Ensure the token contains an "
                "'email' or 'preferred_username' claim."
            )
            log.error(err)
            _sync_summary["errors"].append(err)
            _sync_summary["exit_status"] = "failure"
            sys.exit(1)

        log.info(
            "Resolving account number '%s' via user details for %s",
            account_number, caller_email,
        )
        try:
            data = self._get("/User/Details", params={"email": caller_email})
        except requests.HTTPError as exc:
            err = f"GET /User/Details failed for {caller_email}: {exc}"
            log.error(err)
            _sync_summary["errors"].append(err)
            _sync_summary["exit_status"] = "failure"
            sys.exit(1)

        for acct in (data.get("accounts") or []):
            # Accounts schema: accountId (uuid), accountNumber, accountName
            if (acct.get("accountNumber") == account_number
                    or acct.get("accountName") == account_number):
                uid = acct.get("accountId")
                log.info(
                    "Resolved account '%s' → UUID %s", account_number, uid
                )
                return uid

        err = (
            f"Account '{account_number}' not found in the account list for "
            f"{caller_email}. Ensure the account exists and the PAT owner "
            "is a member of it."
        )
        log.error(err)
        _sync_summary["errors"].append(err)
        _sync_summary["exit_status"] = "failure"
        sys.exit(1)


    # ── AccountUser operations ──────────────────────────────────────────────

    def get_account_users(self, account_uuid: str) -> list[dict]:
        """POST /Account/{AccountId}/Users — all pages."""
        all_users, page = [], 0
        while True:
            data  = self._post(
                f"/Account/{account_uuid}/Users", body=None,
                params={"PageNumber": page, "PageSize": 200},
            )
            items = data.get("users") or []
            all_users.extend(items)
            total = data.get("totalCount", len(all_users))
            if len(all_users) >= total or not items:
                break
            page += 1
        log.info("Account has %d current user(s).", len(all_users))
        return all_users

    def remove_users_from_account(self, account_uuid: str,
                                  user_uuids: list[str]) -> None:
        """
        DELETE /Account/{AccountId}/Users
        Body: array of user UUID strings. Removes multiple users in one call.
        """
        if not user_uuids:
            return
        log.info(
            "  - Removing %d user(s) from account %s",
            len(user_uuids), account_uuid,
        )
        self._delete(f"/Account/{account_uuid}/Users", body=user_uuids)

    def bulk_import_account_users(self, account_uuid: str,
                                  users: list[dict],
                                  remove_others: bool = False) -> dict:
        """PUT /Account/{AccountId}/Users"""
        log.info(
            "Bulk importing %d user(s) into account %s (removeOthers=%s)",
            len(users), account_uuid, remove_others,
        )
        params = {"removeOthers": str(remove_others).lower()}
        if not self.send_email:
            params["SuppressEmail"] = "true"
        if self.dry_run:
            log.info(
                "[DRY-RUN] PUT /Account/%s/Users  params=%s  first 3: %s",
                account_uuid, params, json.dumps(users[:3]),
            )
            return {}
        url  = f"{ANSYS_IAM_BASE}/Account/{account_uuid}/Users"
        resp = requests.put(url, headers=self.headers,
                            json=users, params=params)
        resp.raise_for_status()
        return resp.json() if resp.content else {}

    # ── Group operations ────────────────────────────────────────────────────

    def find_group(self, account_uuid: str,
                   group_name: str) -> dict | None:
        """POST /Account/{AccountId}/Groups?searchValue=..."""
        data  = self._post(
            f"/Account/{account_uuid}/Groups", body=None,
            params={"searchValue": group_name,
                    "PageSize": 50, "PageNumber": 0},
        )
        items = data.get("groups") or []
        for g in items:
            if g.get("name") == group_name:
                log.info(
                    "Found group '%s' → UUID %s", group_name, g.get("id")
                )
                return g
        return None

    def create_group(self, account_uuid: str,
                     group_name: str) -> dict:
        """POST /Account/{AccountId}/Group"""
        log.info(
            "Creating group '%s' under account %s", group_name, account_uuid
        )
        body = {
            "name":        group_name,
            "description": f"Synced from Entra ID group: {group_name}",
            "groupTypeId": 1,   # 1 = General, 0 = Private
        }
        if self.dry_run:
            log.info(
                "[DRY-RUN] POST /Account/%s/Group  body=%s",
                account_uuid, json.dumps(body),
            )
            return {"id": "dry-run-group-uuid"}
        return self._post(f"/Account/{account_uuid}/Group", body=body)

    # ── GroupMember operations ──────────────────────────────────────────────

    def get_group_members(self, account_uuid: str,
                          group_uuid: str) -> list[dict]:
        """POST /Account/{AccountId}/Group/{GroupId}/Members — all pages."""
        all_members, page = [], 0
        while True:
            data  = self._post(
                f"/Account/{account_uuid}/Group/{group_uuid}/Members",
                body=None,
                params={"PageNumber": page, "PageSize": 200},
            )
            items = data.get("members") or []
            all_members.extend(items)
            total = data.get("totalCount", len(all_members))
            if len(all_members) >= total or not items:
                break
            page += 1
        log.info("Group has %d current member(s).", len(all_members))
        return all_members

    def add_members_to_group(self, account_uuid: str, group_uuid: str,
                             member_uuids: list[str]) -> None:
        """
        PUT /Account/{AccountId}/Group/{GroupId}/Members
        Body: MemberInGroupRequestModel { members: [{ memberId: uuid }] }
        """
        if not member_uuids:
            return
        log.info(
            "  + Adding %d member(s) to group %s",
            len(member_uuids), group_uuid,
        )
        body = {"members": [{"memberId": uid} for uid in member_uuids]}
        self._put(
            f"/Account/{account_uuid}/Group/{group_uuid}/Members",
            body=body,
        )

    def remove_members_from_group(self, account_uuid: str, group_uuid: str,
                                  member_uuids: list[str]) -> None:
        """
        DELETE /Account/{AccountId}/Group/{GroupId}/Members
        Body: MemberInGroupRequestModel { members: [{ memberId: uuid }] }
        """
        if not member_uuids:
            return
        log.info(
            "  - Removing %d member(s) from group %s",
            len(member_uuids), group_uuid,
        )
        body = {"members": [{"memberId": uid} for uid in member_uuids]}
        self._delete(
            f"/Account/{account_uuid}/Group/{group_uuid}/Members",
            body=body,
        )

    def import_members_to_group(self, account_uuid: str, group_uuid: str,
                                emails: list[str],
                                remove_others: bool = False) -> dict:
        """PUT /Account/{AccountId}/Group/{GroupId}/Members/Import"""
        log.info(
            "Importing %d user(s) into group %s (removeOthers=%s)",
            len(emails), group_uuid, remove_others,
        )
        body   = {"users": [{"email": e, "isAdmin": False} for e in emails]}
        params = {"removeOthers": str(remove_others).lower()}
        if self.dry_run:
            log.info(
                "[DRY-RUN] PUT /Account/%s/Group/%s/Members/Import  body=%s",
                account_uuid, group_uuid, json.dumps(body),
            )
            return {}
        url  = (
            f"{ANSYS_IAM_BASE}/Account/{account_uuid}"
            f"/Group/{group_uuid}/Members/Import"
        )
        resp = requests.put(url, headers=self.headers,
                            json=body, params=params)
        resp.raise_for_status()
        result = resp.json() if resp.content else {}
        for bucket in ("success", "notModified", "error"):
            entries = result.get(bucket) or []
            if entries:
                log.info("  Group import %s: %d user(s)", bucket, len(entries))
                for e in entries:
                    if e.get("errorMessage"):
                        log.warning("    %s — %s",
                                    e.get("email"), e.get("errorMessage"))
        return result


# ============================================================
# FIELD HELPERS
# ============================================================

def ansys_email(member: dict) -> str | None:
    """
    Extract email from an Ansys API member dict.
    AccountUserDetailsResponse: "email"
    MemberResponseModel:        "email"
    """
    return member.get("email")


def ansys_user_uuid(member: dict) -> str | None:
    """
    Extract UUID from an Ansys API member dict.
    AccountUserDetailsResponse: "userId"
    MemberResponseModel:        "memberId"
    """
    return member.get("userId") or member.get("memberId")


# ============================================================
# SYNC — ACCOUNT
# ============================================================

def sync_to_account(client: AnsysIAMClient, account_uuid: str,
                    entra_members: list[dict]) -> None:
    log.info("=== Syncing to ACCOUNT %s ===", account_uuid)

    # entra_email() already returns lowercased values
    desired_emails: set[str] = {
        e for m in entra_members if (e := entra_email(m))
    }
    log.info("Desired: %d user(s)", len(desired_emails))

    current_members  = client.get_account_users(account_uuid)
    current_by_email: dict[str, str] = {}
    for m in current_members:
        email = ansys_email(m)
        uid   = ansys_user_uuid(m)
        if email and uid:
            current_by_email[email.lower()] = uid

    to_add    = desired_emails - set(current_by_email.keys())
    to_remove = set(current_by_email.keys()) - desired_emails

    # Never remove preserved addresses
    if client.preserve_emails:
        protected = to_remove & client.preserve_emails
        if protected:
            log.info(
                "Skipping removal of %d preserved address(es): %s",
                len(protected), sorted(protected),
            )
        to_remove -= client.preserve_emails

    log.info("To add: %d,  To remove: %d", len(to_add), len(to_remove))
    _sync_summary["added"].extend(sorted(to_add))
    _sync_summary["removed"].extend(sorted(to_remove))

    if to_add:
        # CreateAccountUserRequest: { email, roleId (optional UUID) }
        import_payload = [
            {"email": email}
            for email in to_add
        ]
        client.bulk_import_account_users(
            account_uuid, import_payload, remove_others=False
        )

    if to_remove:
        uuids_to_remove = [current_by_email[e] for e in to_remove]
        client.remove_users_from_account(account_uuid, uuids_to_remove)

    log.info("Account sync complete.")


# ============================================================
# SYNC — GROUP
# ============================================================

def sync_to_group(client: AnsysIAMClient, account_uuid: str,
                  group_name: str,
                  entra_members: list[dict]) -> None:
    log.info(
        "=== Syncing to GROUP '%s' under account %s ===",
        group_name, account_uuid,
    )

    group = client.find_group(account_uuid, group_name)
    if group is None:
        log.info("Group '%s' not found — creating it.", group_name)
        group = client.create_group(account_uuid, group_name)

    # find_group returns GroupDetailsResponse (field: "id")
    # create_group returns CreateGroupResponse (field: "groupId")
    group_uuid = group.get("groupId") or group.get("id")
    if not group_uuid or group_uuid == "dry-run-group-uuid":
        if client.dry_run:
            log.info("[DRY-RUN] Skipping member sync — no real group UUID.")
            return
        err = f"Could not resolve group UUID. Response: {group}"
        log.error(err)
        _sync_summary["errors"].append(err)
        _sync_summary["exit_status"] = "failure"
        sys.exit(1)
    log.info("Using group UUID: %s", group_uuid)

    desired_emails: set[str] = {
        e for m in entra_members if (e := entra_email(m))
    }
    log.info("Desired: %d member(s)", len(desired_emails))

    current_members  = client.get_group_members(account_uuid, group_uuid)
    current_by_email: dict[str, str] = {}
    for m in current_members:
        if m.get("isGroupAsMember"):
            continue
        email = ansys_email(m)
        uid   = ansys_user_uuid(m)
        if email and uid:
            current_by_email[email.lower()] = uid

    to_add    = desired_emails - set(current_by_email.keys())
    to_remove = set(current_by_email.keys()) - desired_emails

    # Never remove preserved addresses
    if client.preserve_emails:
        protected = to_remove & client.preserve_emails
        if protected:
            log.info(
                "Skipping removal of %d preserved address(es): %s",
                len(protected), sorted(protected),
            )
        to_remove -= client.preserve_emails

    log.info("To add: %d,  To remove: %d", len(to_add), len(to_remove))
    _sync_summary["added"].extend(sorted(to_add))
    _sync_summary["removed"].extend(sorted(to_remove))

    if to_add:
        client.import_members_to_group(
            account_uuid, group_uuid, list(to_add), remove_others=False
        )

    if to_remove:
        uuids_to_remove = [current_by_email[e] for e in to_remove]
        client.remove_members_from_group(
            account_uuid, group_uuid, uuids_to_remove
        )

    log.info("Group sync complete.")


# ============================================================
# CLI
# ============================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Sync an Entra ID security group → "
            "Ansys ID Portal account or group."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Core sync arguments
    g = parser.add_argument_group("Sync")
    g.add_argument("--entra-domain",   required=True,
                   help="Entra ID tenant domain, e.g. contoso.onmicrosoft.com")
    g.add_argument("--entra-group",    required=True,
                   help="Display name of the Entra ID security group")
    g.add_argument("--target-type",    required=True,
                   choices=["account", "group"])
    g.add_argument("--account-number", required=True,
                   help="Ansys ID Portal account number (must already exist)")
    g.add_argument("--group-name",     default=None,
                   help="Ansys group name (required when --target-type is group)")
    g.add_argument("--dry-run",        action="store_true",
                   help="Print planned changes without applying them")
    g.add_argument("--send-email", action="store_true",
                   help=(
                       "Send Ansys ID invitation emails when adding users. "
                       "Emails are suppressed by default."
                   ))
    g.add_argument(
        "--reset-credentials", action="store_true",
        help=(
            "Clear all cached credentials (Entra client ID, client secret, "
            "and Ansys PAT) from the OS keychain and prompt for new values."
        ),
    )
    g.add_argument(
        "--preserve-email", default=None,
        help=(
            "Comma-separated email address(es) to never remove from the "
            "Ansys target, even if absent from the Entra group. Useful for "
            "protecting the PAT owner or service accounts from being synced "
            "out. Example: michael.long@ansys.com,svc-account@ansys.com"
        ),
    )

    # Ansys B2C credentials
    b = parser.add_argument_group("Ansys B2C Credentials")
    b.add_argument(
        "--ansys-pat", default=None,
        help=(
            "Ansys ID Portal Personal Access Token used for B2C authentication. "
            "If omitted, the script will prompt interactively. "
            "Avoid passing on the command line in production -- "
            "use the interactive prompt instead."
        ),
    )

    # Entra credentials
    c = parser.add_argument_group(
        "Entra ID Credentials",
        "If omitted, the script will prompt interactively. "
        "The client secret prompt does not echo to the screen.",
    )
    c.add_argument(
        "--entra-client-id", default=None,
        help="App registration client ID for Graph API access",
    )
    c.add_argument(
        "--entra-client-secret", default=None,
        help=(
            "App registration client secret for Graph API access. "
            "Avoid passing on the command line in production — "
            "use the interactive prompt instead."
        ),
    )

    # Logging / notification
    n = parser.add_argument_group("Logging & Notification")
    n.add_argument("--log-file",       default=None,
                   help="Path to append the full run log to")
    n.add_argument("--smtp-relay",     default=None,
                   help="Hostname or IP of the SMTP relay/gateway")
    n.add_argument("--smtp-port",      type=int, default=25,
                   help="SMTP port (default: 25)")
    n.add_argument("--smtp-from",      default=None,
                   help="Envelope From address for all outbound emails")
    n.add_argument("--log-recipient",  default=None,
                   help="Comma-separated address(es) for the full log email")
    n.add_argument("--warn-recipient", default=None,
                   help="Comma-separated address(es) for warnings/errors email")

    return parser.parse_args()


def _split_recipients(value: str | None) -> list[str]:
    if not value:
        return []
    return [r.strip() for r in value.split(",") if r.strip()]


def _parse_preserve_emails(value: str | None) -> set[str]:
    """Parse --preserve-email into a lowercase set for O(1) lookup."""
    if not value:
        return set()
    return {e.strip().lower() for e in value.split(",") if e.strip()}


async def main() -> None:
    args       = parse_args()
    log_buffer = setup_logging(args.log_file)
    configure_environment()

    _sync_summary["started_at"]  = datetime.now(timezone.utc).isoformat()
    _sync_summary["entra_group"] = args.entra_group
    _sync_summary["target_type"] = args.target_type
    _sync_summary["target_name"] = args.group_name or args.account_number
    _sync_summary["dry_run"]     = args.dry_run

    log_recipients  = _split_recipients(args.log_recipient)
    warn_recipients = _split_recipients(args.warn_recipient)

    if (log_recipients or warn_recipients) and not args.smtp_relay:
        log.warning(
            "Email recipients specified but --smtp-relay not provided. "
            "No emails will be sent."
        )
    if (log_recipients or warn_recipients) and not args.smtp_from:
        log.warning(
            "Email recipients specified but --smtp-from not provided. "
            "No emails will be sent."
        )

    if args.target_type == "group" and not args.group_name:
        err = "--group-name is required when --target-type is 'group'."
        log.error(err)
        _sync_summary["errors"].append(err)
        _sync_summary["exit_status"] = "failure"
        _sync_summary["finished_at"] = datetime.now(timezone.utc).isoformat()
        send_warn_email(args.smtp_relay, args.smtp_port,
                        args.smtp_from, warn_recipients)
        send_log_email(args.smtp_relay, args.smtp_port, args.smtp_from,
                       log_recipients, log_buffer)
        sys.exit(1)

    try:
        # Clear cached credentials if requested
        if args.reset_credentials:
            reset_cached_credentials()

        # Resolve Entra credentials (prompt if not on command line)
        client_id, client_secret = resolve_entra_credentials(args)

        # 1. Fetch Entra members
        entra_members = get_entra_group_members(
            args.entra_domain, args.entra_group,
            client_id, client_secret,
        )
        if not entra_members:
            log.warning(
                "Entra group returned 0 eligible members. Nothing to sync."
            )
            _sync_summary["exit_status"] = "no-op"
            return

        # 2. Obtain Ansys token
        ansys_token, caller_email = await get_ansys_access_token(
            CLIENT_ID, ANSYS_SCOPE,
            pat_arg=args.ansys_pat,
        )
        if caller_email:
            log.info("Authenticated as: %s", caller_email)
        else:
            log.warning(
                "Could not extract email from token — "
                "account UUID resolution will fail."
            )
        client = AnsysIAMClient(
            ansys_token,
            dry_run=args.dry_run,
            send_email=args.send_email,
            preserve_emails=_parse_preserve_emails(args.preserve_email),
        )

        # 3. Resolve account number → UUID via authenticated user's account list
        account_uuid = client.find_account_uuid(
            args.account_number, caller_email=caller_email
        )

        # 4. Sync
        if args.target_type == "account":
            sync_to_account(client, account_uuid, entra_members)
        else:
            sync_to_group(
                client, account_uuid, args.group_name, entra_members
            )

    except Exception as exc:
        err = f"Unhandled exception: {exc}"
        log.error(err, exc_info=True)
        _sync_summary["errors"].append(err)
        _sync_summary["exit_status"] = "failure"

    finally:
        _sync_summary["finished_at"] = datetime.now(timezone.utc).isoformat()
        log.info(
            "Sync finished — status: %s  added: %d  removed: %d  "
            "skipped: %d  errors: %d",
            _sync_summary["exit_status"],
            len(_sync_summary["added"]),
            len(_sync_summary["removed"]),
            len(_sync_summary["skipped"]),
            len(_sync_summary["errors"]),
        )
        send_warn_email(
            args.smtp_relay, args.smtp_port,
            args.smtp_from, warn_recipients,
        )
        send_log_email(
            args.smtp_relay, args.smtp_port,
            args.smtp_from, log_recipients,
            log_buffer,
        )


if __name__ == "__main__":
    asyncio.run(main())