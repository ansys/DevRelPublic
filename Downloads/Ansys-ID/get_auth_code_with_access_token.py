"""
Dynamic Azure AD B2C Token Management + API Caller
---------------------------------------------------
Flow:
  1️ Try cached access token.
  2 If expired, refresh using refresh token.
  3 If no refresh token, ask for PAT (and optional client IP).
  4 Call target API with valid bearer token.

Dependencies:
  pip install msal httpx requests

Run
  py get_auth_code_with_access_token.py
"""

import sys
import json
import requests
import httpx
import msal
from pathlib import Path
from urllib.parse import urlparse, parse_qs

# ============================================================
# CONFIGURATION
# ============================================================
AADB2C_ENV = "ansysaccount"
TENANT_NAME = f"{AADB2C_ENV}.onmicrosoft.com"
POLICY = "B2C_1A_ANSYSID_SIGNUP_SIGNIN"

CLIENT_ID = "28982bbf-f354-4e48-8bfb-e542d44c588c"
REDIRECT_URI = "https://login.microsoftonline.com/common/oauth2/nativeclient"
API_SCOPE = f"https://{TENANT_NAME}/LicensingAsAService/LicensingAsAService openid offline_access"

TOKEN_CACHE_FILE = Path(".token_cache_prod.json")
API_URL = "https://example.ansys.com/v1/products?verbose=true&accountId=test_api_feb_02_0&all=true"

# ============================================================
# TOKEN CACHE HELPERS
# ============================================================
def load_token_cache():
    cache = msal.SerializableTokenCache()
    if TOKEN_CACHE_FILE.exists():
        cache.deserialize(TOKEN_CACHE_FILE.read_text())
    return cache


def save_token_cache(cache):
    if cache.has_state_changed:
        TOKEN_CACHE_FILE.write_text(cache.serialize())


# ============================================================
# AUTH HELPERS
# ============================================================
def build_authorization_url(client_id: str, pat: str, client_ip: str = None) -> str:
    """Construct the Azure AD B2C authorization URL using PAT as id_token_hint."""
    base = (
        f"https://{AADB2C_ENV}.b2clogin.com/"
        f"{TENANT_NAME}/oauth2/v2.0/authorize"
        f"?p={POLICY}"
        f"&client_id={client_id}"
        f"&nonce=defaultNonce"
        f"&redirect_uri={REDIRECT_URI}"
        f"&scope=openid"
        f"&response_type=code"
        f"&id_token_hint={pat}"
    )
    if client_ip:
        base += f"&client_ip={client_ip}"
    return base


async def get_auth_code(client: httpx.AsyncClient, auth_url: str) -> str:
    """Obtain authorization code automatically from Azure AD B2C redirect."""
    resp = await client.get(auth_url, follow_redirects=False)
    location = resp.headers.get("Location")
    if not location:
        raise RuntimeError("No redirect received — check PAT validity.")
    parsed = urlparse(location)
    params = parse_qs(parsed.query)
    if "code" in params:
        return params["code"][0]
    if "error" in params:
        raise RuntimeError(params.get("error_description", ["Unknown error"])[0])
    raise RuntimeError("Authorization code not received.")


# ============================================================
# TOKEN ACQUISITION
# ============================================================
def acquire_token_silent(app):
    accounts = app.get_accounts()
    if accounts:
        result = app.acquire_token_silent([API_SCOPE], account=accounts[0])
        if result and "access_token" in result:
            # print("✅ Access token retrieved silently from cache.")
            return result
    return None


def exchange_code_for_tokens(app, code):
    return app.acquire_token_by_authorization_code(
        code=code,
        scopes=[API_SCOPE],
        redirect_uri=REDIRECT_URI
    )


# ============================================================
# MAIN FUNCTION
# ============================================================
async def get_access_token():
    cache = load_token_cache()
    authority = f"https://{AADB2C_ENV}.b2clogin.com/{TENANT_NAME}/{POLICY}"
    app = msal.PublicClientApplication(client_id=CLIENT_ID, authority=authority, token_cache=cache)

    # Step 1: Try from cache
    result = acquire_token_silent(app)
    if result:
        save_token_cache(cache)
        return result["access_token"]

    # Step 2: Ask for PAT
    print("⚠️  No valid token found. Please provide your PAT for authentication.")
    pat = input("Enter your PAT: ").strip()
    client_ip = input("Optional - Enter your Client IP (or leave blank): ").strip() or None

    print("🔐 Authenticating using PAT...")
    auth_url = build_authorization_url(CLIENT_ID, pat, client_ip)
    print(auth_url)
    async with httpx.AsyncClient(timeout=15) as client:
        try:
            code = await get_auth_code(client, auth_url)
            print("✅ Authorization code retrieved successfully.")
        except Exception as e:
            print(f"❌ Failed to obtain authorization code: {e}")
            sys.exit(1)

    # Step 3: Exchange for tokens
    result = exchange_code_for_tokens(app, code)
    if "access_token" in result:
        print("✅ Access token acquired and cached.")
        save_token_cache(cache)
        return result["access_token"]
    else:
        print("❌ Authentication failed:", result.get("error_description"))
        sys.exit(1)


# ============================================================
# CALL API
# ============================================================
async def call_api():
    print(f"\n🚀 Calling API: {API_URL}\n")
    response = requests.get(API_URL, headers={"Authorization": f"Bearer {await get_access_token()}"})
    print("Status:", response.status_code)
    try:
        print(json.dumps(response.json(), indent=2))
    except Exception:
        print(response.text)


# ============================================================
if __name__ == "__main__":
    import asyncio
    asyncio.run(call_api())
