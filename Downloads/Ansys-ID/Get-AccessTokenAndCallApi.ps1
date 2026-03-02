<#
Dynamic Azure AD B2C Token Management + API Caller (PowerShell)
---------------------------------------------------------------

✔ Uses correct Azure AD B2C token endpoint
✔ Proper Authorization Code → Token exchange
✔ Uses PAT as id_token_hint
✔ Saves & loads refresh tokens
✔ Calls downstream API with bearer token

Install-Module MSAL.PS -Scope CurrentUser

#>

# ==========================
# CONFIG
# ==========================
$AADB2CEnv = "ansysaccount"
$Tenant = "$AADB2CEnv.onmicrosoft.com"
$Policy = "B2C_1A_ANSYSID_SIGNUP_SIGNIN"

$ClientId = "28982bbf-f354-4e48-8bfb-e542d44c588c"
$RedirectUri = "https://login.microsoftonline.com/common/oauth2/nativeclient"
$ApiScope = "https://$Tenant/LicensingAsAService/LicensingAsAService offline_access openid"

$ApiUrl = "https://example.ansys.com/v1/products?verbose=true&accountId=test_api_feb_02_0&all=true"
$TokenCacheFile = ".\ps_token_cache.json"


# ==========================
# Load/Save Token Cache
# ==========================
function Load-TokenCache {
    if (Test-Path $TokenCacheFile) {
        return Get-Content $TokenCacheFile -Raw | ConvertFrom-Json
    }
    return $null
}

function Save-TokenCache($token) {
    $token | ConvertTo-Json -Depth 10 | Out-File $TokenCacheFile
}


# ==========================
# Silent Token Retrieval
# ==========================
function Try-SilentToken {
    $cache = Load-TokenCache
    if ($cache -and $cache.refresh_token) {
        Write-Host "Refreshing access token using cached refresh token..."

        $tokenUrl = "https://$AADB2CEnv.b2clogin.com/$Tenant/oauth2/v2.0/token?p=$Policy"

        $body = @{
            client_id    = $ClientId
            grant_type   = "refresh_token"
            refresh_token = $cache.refresh_token
            scope        = $ApiScope
            redirect_uri = $RedirectUri
        }

        try {
            $response = Invoke-RestMethod -Uri $tokenUrl -Method POST -Body $body -ContentType "application/x-www-form-urlencoded"
            Save-TokenCache $response
            Write-Host "Token refreshed successfully."
            return $response.access_token
        }
        catch {
            Write-Host "Failed to refresh token: $($_.Exception.Message)"
        }
    }
    return $null
}


# ==========================
# Build Authorization URL
# ==========================
function Build-AuthUrl($pat, $clientIp) {

    $url = "https://$AADB2CEnv.b2clogin.com/$Tenant/oauth2/v2.0/authorize" +
        "?p=$Policy" +
        "&client_id=$ClientId" +
        "&nonce=defaultNonce" +
        "&redirect_uri=$RedirectUri" +
        "&scope=openid" +
        "&response_type=code" +
        "&id_token_hint=$pat"

    if ($clientIp) { $url += "&client_ip=$clientIp" }

    return $url
}


# ==========================
# Extract Authorization Code
# ==========================
function Get-AuthCode($authUrl) {

    Write-Host "Requesting authorization code..."

    $resp = Invoke-WebRequest -Uri $authUrl -MaximumRedirection 0 -ErrorAction SilentlyContinue
    $location = $resp.Headers["Location"]

    if (-not $location) {
        throw "Failed: No redirect received (check PAT validity)."
    }

    Write-Host "Redirect received."

    $uri = [System.Uri]$location
    $qs = [System.Web.HttpUtility]::ParseQueryString($uri.Query)

    if ($qs["code"]) {
        Write-Host "Authorization code retrieved."
        return $qs["code"]
    }

    if ($qs["error_description"]) {
        throw "Azure AD B2C Error: $($qs['error_description'])"
    }

    throw "Authorization code missing in redirect response."
}


# ==========================
# Exchange Code → Token (Corrected B2C endpoint)
# ==========================
function Exchange-CodeForToken($code) {

    Write-Host "Exchanging authorization code for access token..."

    $tokenUrl = "https://$AADB2CEnv.b2clogin.com/$Tenant/oauth2/v2.0/token?p=$Policy"

    $body = @{
        grant_type   = "authorization_code"
        code         = $code
        client_id    = $ClientId
        redirect_uri = $RedirectUri
        scope        = $ApiScope
    }

    try {
        $response = Invoke-RestMethod -Uri $tokenUrl -Method POST -Body $body -ContentType "application/x-www-form-urlencoded"
        Save-TokenCache $response
        Write-Host "Access token acquired & cached."
        return $response.access_token
    }
    catch {
        Write-Host "Error exchanging code for token: $($_.Exception.Message)"
        throw
    }
}


# ==========================
# MAIN TOKEN FUNCTION
# ==========================
function Get-AccessToken {

    # Step 1 — Try silent refresh
    $silentToken = Try-SilentToken
    if ($silentToken) { return $silentToken }

    # Step 2 — Get PAT
    Write-Host "No cached token available. PAT is required."
    $pat = Read-Host "Enter your PAT"
    $clientIp = Read-Host "Optional - Enter your Client IP (or press Enter)"

    # Step 3 — Build URL and get code
    $authUrl = Build-AuthUrl $pat $clientIp
    Write-Host "`nOpen this URL for debugging if needed:"
    Write-Host $authUrl

    $code = Get-AuthCode $authUrl

    # Step 4 — Exchange code
    return Exchange-CodeForToken $code
}


# ==========================
# CALL TARGET API
# ==========================
function Call-Api {
    $token = Get-AccessToken

    Write-Host "`nCalling API..."
    $headers = @{ Authorization = "Bearer $token" }

    try {
        $resp = Invoke-RestMethod -Uri $ApiUrl -Headers $headers -Method GET
        Write-Host "Status: Success"
        $resp | ConvertTo-Json -Depth 10
    }
    catch {
        Write-Host "API Error:"
        Write-Host $_.Exception.Message
    }
}


# ==========================
# RUN MAIN
# ==========================
Call-Api
