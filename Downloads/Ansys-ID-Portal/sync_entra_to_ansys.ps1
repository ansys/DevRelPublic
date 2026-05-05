<#
.SYNOPSIS
    Synchronizes members of an Azure Entra ID security group to an Ansys ID
    Portal account or group.

.DESCRIPTION
    Entra ID -> Ansys ID Portal Sync
    ---------------------------------
    Reads the membership of a named Entra ID security group via the Microsoft
    Graph API, then brings the target Ansys ID Portal account or group into
    exact alignment -- adding missing users and removing extras.

    - Accounts must already exist; the script resolves the account number to
      a UUID via the authenticated user's own account list.
    - Groups are created automatically if they do not exist.
    - Members without a primary email (mail) attribute in Entra are skipped
      with a warning.
    - Preserved email addresses are never removed from the target.

    No modules beyond the PowerShell standard library are required.
    All HTTP calls use Invoke-RestMethod / Invoke-WebRequest.

    API Reference: https://iam.ansys.com/swagger/AnsysId/swagger.json

.PARAMETER EntraDomain
    Entra ID tenant domain, e.g. contoso.onmicrosoft.com

.PARAMETER EntraGroup
    Display name of the Entra ID security group to sync from.

.PARAMETER EntraClientId
    App registration client ID for Microsoft Graph API access.
    Prompted interactively if not supplied.

.PARAMETER EntraClientSecret
    App registration client secret for Microsoft Graph API access.
    Prompted interactively (hidden) if not supplied.
    Avoid passing on the command line in production.

.PARAMETER TargetType
    Either 'account' or 'group'.

.PARAMETER AccountNumber
    Ansys ID Portal account number (must already exist).

.PARAMETER GroupName
    Ansys ID Portal group name. Required when TargetType is 'group'.
    Created automatically if it does not exist.


.PARAMETER AnsysPat
    Ansys ID Portal Personal Access Token for B2C authentication.
    Prompted interactively (hidden) if not supplied.
    Avoid passing on the command line in production.


.PARAMETER DryRun
    Switch. Log planned changes without writing anything to Ansys ID Portal.

.PARAMETER SendEmail
    Switch. Send Ansys ID invitation emails when adding users.
    Emails are suppressed by default.

.PARAMETER PreserveEmail
    Comma-separated email address(es) never removed from the target, even if
    absent from the Entra group. Protects the PAT owner or service accounts.

.PARAMETER LogFile
    Path to a file where the full run log is appended. Directory must exist.

.PARAMETER SmtpRelay
    Hostname or IP of the SMTP relay for notification emails.

.PARAMETER SmtpPort
    SMTP port. Default: 25.

.PARAMETER SmtpFrom
    Envelope From address for all outbound emails.

.PARAMETER LogRecipient
    Comma-separated address(es) for the full log email after every run.

.PARAMETER WarnRecipient
    Comma-separated address(es) for a warnings/errors summary email.
    Only sent when there is something actionable to report.

.EXAMPLE
    .\sync_entra_to_ansys.ps1 `
        -EntraDomain   'contoso.onmicrosoft.com' `
        -EntraGroup    'Ansys License Users' `
        -TargetType    account `
        -AccountNumber 'MDLTestAccount' `

.EXAMPLE
    .\sync_entra_to_ansys.ps1 `
        -EntraDomain       'contoso.onmicrosoft.com' `
        -EntraGroup        'Ansys License Users' `
        -EntraClientId     '2874d43f-518b-4271-b1ef-04527842b3f4' `
        -EntraClientSecret 'your-secret' `
        -TargetType        group `
        -AccountNumber     'MDLTestAccount' `
        -GroupName         'AnsysLicenseGroup' `
        -PreserveEmail     'admin@contoso.com' `
        -LogFile           'C:\Logs\ansys_sync.log' `
        -SmtpRelay         'smtp.contoso.com' `
        -SmtpFrom          'ansys-sync@contoso.com' `
        -LogRecipient      'it-admin@contoso.com' `
        -WarnRecipient     'it-alerts@contoso.com'
#>

[CmdletBinding()]
param (
    # Required sync parameters
    [Parameter(Mandatory)][string] $EntraDomain,
    [Parameter(Mandatory)][string] $EntraGroup,
    [Parameter(Mandatory)][ValidateSet('account','group')][string] $TargetType,
    [Parameter(Mandatory)][string] $AccountNumber,

    # Optional sync parameters
    [string] $GroupName        = '',
    [switch] $DryRun,
    [switch] $SendEmail,
    [string] $PreserveEmail    = '',

    # Entra ID credentials (prompted if omitted)
    [string] $EntraClientId     = '',
    [string] $EntraClientSecret = '',

    # Ansys B2C credentials
    [string] $AnsysPat      = '',

    # Logging and notification
    [string] $LogFile       = '',
    [string] $SmtpRelay     = '',
    [int]    $SmtpPort      = 25,
    [string] $SmtpFrom      = '',
    [string] $LogRecipient  = '',
    [string] $WarnRecipient = ''
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

# ============================================================
# LOGGING
# ============================================================

$script:LogLines      = [System.Collections.Generic.List[string]]::new()
$script:WarnLines     = [System.Collections.Generic.List[string]]::new()
$script:LogFileHandle = $null

function Write-Log {
    param(
        [ValidateSet('INFO','WARNING','ERROR','DEBUG')][string] $Level = 'INFO',
        [string] $Message
    )
    $ts   = (Get-Date).ToUniversalTime().ToString('yyyy-MM-ddTHH:mm:ss')
    $line = '{0} [{1}] {2}' -f $ts, $Level.PadRight(7), $Message

    $script:LogLines.Add($line)

    if ($Level -eq 'WARNING' -or $Level -eq 'ERROR') {
        $script:WarnLines.Add($line)
    }

    if ($script:LogFileHandle) {
        $script:LogFileHandle.WriteLine($line)
        $script:LogFileHandle.Flush()
    }

    switch ($Level) {
        'ERROR'   { Write-Host $line -ForegroundColor Red    }
        'WARNING' { Write-Host $line -ForegroundColor Yellow }
        'DEBUG'   { <# suppress debug from console #>        }
        default   { Write-Host $line                         }
    }
}

function Open-LogFile {
    param([string] $Path)
    if (-not $Path) { return }
    $script:LogFileHandle = [System.IO.StreamWriter]::new(
        $Path, $true, [System.Text.Encoding]::UTF8
    )
    Write-Log INFO "Log file: $Path"
}

function Close-LogFile {
    if ($script:LogFileHandle) {
        $script:LogFileHandle.Close()
        $script:LogFileHandle = $null
    }
}

# ============================================================
# SYNC SUMMARY
# ============================================================

$script:Summary = @{
    StartedAt  = $null
    FinishedAt = $null
    EntraGroup = $EntraGroup
    TargetType = $TargetType
    TargetName = if ($GroupName) { $GroupName } else { $AccountNumber }
    Added      = [System.Collections.Generic.List[string]]::new()
    Removed    = [System.Collections.Generic.List[string]]::new()
    Skipped    = [System.Collections.Generic.List[string]]::new()
    Errors     = [System.Collections.Generic.List[string]]::new()
    DryRun     = $DryRun.IsPresent
    ExitStatus = 'success'
}

# ============================================================
# SMTP HELPERS
# ============================================================

function Split-Recipients {
    param([string] $Value)
    if (-not $Value) { return @() }
    return @($Value -split ',' | ForEach-Object { $_.Trim() } | Where-Object { $_ })
}

function Send-SyncEmail {
    param(
        [string[]] $To,
        [string]   $Subject,
        [string]   $Body
    )
    if (-not $SmtpRelay -or -not $To -or -not $SmtpFrom) { return }
    try {
        $msg            = [System.Net.Mail.MailMessage]::new()
        $msg.From       = $SmtpFrom
        $msg.Subject    = $Subject
        $msg.Body       = $Body
        $msg.IsBodyHtml = $false
        foreach ($addr in $To) { $msg.To.Add($addr) }
        $smtp = [System.Net.Mail.SmtpClient]::new($SmtpRelay, $SmtpPort)
        $smtp.Send($msg)
        $msg.Dispose()
        $smtp.Dispose()
        Write-Log INFO "Email sent to $($To -join ', ') -- subject: $Subject"
    }
    catch {
        Write-Log WARNING "Failed to send email to $($To -join ', '): $_"
    }
}

function Send-LogEmail {
    param([string[]] $Recipients)
    if (-not $Recipients -or $Recipients.Count -eq 0) { return }

    $s      = $script:Summary
    $tag    = if ($s.DryRun) { '[DRY-RUN] ' } else { '' }
    $ts     = (Get-Date).ToUniversalTime().ToString('yyyy-MM-dd HH:mm UTC')
    $status = $s.ExitStatus.ToUpper()
    $dr     = if ($s.DryRun) { ' (DRY-RUN)' } else { '' }
    $subj   = '{0}Ansys Sync {1} -- {2} -> {3}:{4} -- {5}' -f $tag, $status, $s.EntraGroup, $s.TargetType, $s.TargetName, $ts

    $sep    = '=' * 60
    $added  = ($s.Added   | Sort-Object) -join ', '
    $removed= ($s.Removed | Sort-Object) -join ', '
    $skipped= ($s.Skipped | Sort-Object) -join ', '
    $log    = $script:LogLines -join "`n"

    $body = "Ansys ID Portal Sync -- Full Log`n"
    $body += "$sep`n"
    $body += "Started :  $($s.StartedAt)`n"
    $body += "Finished:  $($s.FinishedAt)`n"
    $body += "Status  :  $status$dr`n"
    $body += "Source  :  Entra group '$($s.EntraGroup)'`n"
    $body += "Target  :  $($s.TargetType) '$($s.TargetName)'`n"
    $body += "`n"
    $body += "Added   : $($s.Added.Count)  $added`n"
    $body += "Removed : $($s.Removed.Count)  $removed`n"
    $body += "Skipped : $($s.Skipped.Count)  $skipped`n"
    $body += "Errors  : $($s.Errors.Count)`n"
    $body += "`n"
    $body += "$sep`n"
    $body += "FULL LOG`n"
    $body += "$sep`n"
    $body += $log

    Send-SyncEmail -To $Recipients -Subject $subj -Body $body
}

function Send-WarnEmail {
    param([string[]] $Recipients)
    if (-not $Recipients -or $Recipients.Count -eq 0) { return }
    if ($script:WarnLines.Count -eq 0 -and $script:Summary.Errors.Count -eq 0) {
        Write-Log INFO "No warnings or errors -- skipping alert email."
        return
    }

    $s    = $script:Summary
    $tag  = if ($s.DryRun) { '[DRY-RUN] ' } else { '' }
    $ts   = (Get-Date).ToUniversalTime().ToString('yyyy-MM-dd HH:mm UTC')
    $dr   = if ($s.DryRun) { ' (DRY-RUN)' } else { '' }
    $subj = '{0}[WARN] Ansys Sync WARNINGS -- {1} -> {2}:{3} -- {4}' -f $tag, $s.EntraGroup, $s.TargetType, $s.TargetName, $ts

    $sep  = '=' * 60
    $body = "Ansys ID Portal Sync -- Warnings & Errors Summary`n"
    $body += "$sep`n"
    $body += "Run time : $($s.StartedAt) -> $($s.FinishedAt)`n"
    $body += "Status   : $($s.ExitStatus.ToUpper())$dr`n"
    $body += "Source   : Entra group '$($s.EntraGroup)'`n"
    $body += "Target   : $($s.TargetType) '$($s.TargetName)'`n"
    $body += "`n"

    if ($s.Skipped.Count -gt 0) {
        $body += "SKIPPED USERS ($($s.Skipped.Count)) -- no primary email (mail) attribute in Entra:`n"
        foreach ($entry in ($s.Skipped | Sort-Object)) { $body += "  * $entry`n" }
        $body += "`n"
    }
    if ($s.Errors.Count -gt 0) {
        $body += "ERRORS ($($s.Errors.Count)):`n"
        foreach ($err in $s.Errors) { $body += "  * $err`n" }
        $body += "`n"
    }
    if ($script:WarnLines.Count -gt 0) {
        $body += "LOG WARNINGS/ERRORS ($($script:WarnLines.Count)):`n"
        foreach ($wl in $script:WarnLines) { $body += "  $wl`n" }
    }

    Send-SyncEmail -To $Recipients -Subject $subj -Body $body
}

# ============================================================
# ENVIRONMENT CONFIGURATION
# ============================================================

# Production endpoints -- hardcoded.
$B2cEnv         = 'ansysaccount'
$TenantName     = 'ansysaccount.onmicrosoft.com'
$AnsysIamBase   = 'https://iam.ansys.com'
$TokenCacheFile = '.token_cache.json'
$ResolvedScope  = 'https://ansysaccount.onmicrosoft.com/AnsysID/Authentication'
$Policy         = 'B2C_1A_ANSYSID_SIGNUP_SIGNIN'
$RedirectUri    = 'https://login.microsoftonline.com/common/oauth2/nativeclient'
$GraphBase      = 'https://graph.microsoft.com/v1.0'

Write-Log INFO "B2C: https://$B2cEnv.b2clogin.com  ID Portal: $AnsysIamBase  Scope: $ResolvedScope"

# ============================================================
# CREDENTIAL HELPERS
# ============================================================

function Resolve-Credential {
    param(
        [string] $Value,
        [string] $Prompt,
        [string] $Label,
        [switch] $Secret
    )
    if ($Value) { return $Value }
    if ($Secret) {
        $sec  = Read-Host -Prompt $Prompt -AsSecureString
        $bstr = [System.Runtime.InteropServices.Marshal]::SecureStringToBSTR($sec)
        try {
            $plain = [System.Runtime.InteropServices.Marshal]::PtrToStringAuto($bstr)
        }
        finally {
            [System.Runtime.InteropServices.Marshal]::ZeroFreeBSTR($bstr)
        }
        if (-not $plain) { throw "$Label is required." }
        return $plain
    }
    $val = (Read-Host -Prompt $Prompt).Trim()
    if (-not $val) { throw "$Label is required." }
    return $val
}

# ============================================================
# SAFE PROPERTY HELPER
# ============================================================

function Get-SafeProperty {
    <#
    .SYNOPSIS
        Safely retrieves a property from a PSObject without throwing under
        Set-StrictMode -Version Latest when the property does not exist.
    .PARAMETER Obj
        The object to read from.
    .PARAMETER Name
        The property name.
    .PARAMETER Default
        Value to return if the property does not exist. Default: $null.
    #>
    param([object] $Obj, [string] $Name, [object] $Default = $null)
    if ($null -eq $Obj) { return $Default }
    $prop = $Obj.PSObject.Properties[$Name]
    if ($prop) { return $prop.Value }
    return $Default
}

# ============================================================
# ENTRA ID -- GRAPH API
# ============================================================

function Get-GraphToken {
    param([string] $ClientId, [string] $ClientSecret)
    Write-Log INFO "Obtaining Microsoft Graph API token..."
    $tokenUrl = "https://login.microsoftonline.com/$EntraDomain/oauth2/v2.0/token"
    $body = @{
        grant_type    = 'client_credentials'
        client_id     = $ClientId
        client_secret = $ClientSecret
        scope         = 'https://graph.microsoft.com/.default'
    }
    try {
        $resp = Invoke-RestMethod -Uri $tokenUrl -Method POST -Body $body
        return $resp.access_token
    }
    catch {
        throw "Failed to obtain Graph API token: $_"
    }
}

function Get-EntraGroupMembers {
    param([string] $GraphToken)

    $headers = @{ Authorization = "Bearer $GraphToken" }

    Write-Log INFO "Fetching Entra ID group '$EntraGroup' in domain '$EntraDomain'"

    # Find group by display name
    $filterStr = "displayName eq '$EntraGroup' and securityEnabled eq true"
    $encoded   = [Uri]::EscapeDataString($filterStr)
    $uri       = "$GraphBase/groups?`$filter=$encoded&`$select=id,displayName"
    $resp      = Invoke-RestMethod -Uri $uri -Headers $headers -Method GET
    $groups    = $resp.value

    if (-not $groups -or $groups.Count -eq 0) {
        throw "Entra ID security group '$EntraGroup' not found."
    }
    if ($groups.Count -gt 1) {
        Write-Log WARNING "Multiple groups matched '$EntraGroup'; using first result."
    }
    $groupId = $groups[0].id
    Write-Log INFO "Found Entra group id: $groupId"

    # Page through transitive members
    $members = [System.Collections.Generic.List[object]]::new()
    $nextUri = "$GraphBase/groups/$groupId/transitiveMembers?`$select=displayName,mail,userPrincipalName,id&`$top=999"

    while ($nextUri) {
        $page = Invoke-RestMethod -Uri $nextUri -Headers $headers -Method GET
        $pageValue = Get-SafeProperty $page 'value' @()
        foreach ($m in $pageValue) {
            if ($m.'@odata.type' -eq '#microsoft.graph.user') {
                $members.Add($m)
            }
        }
        $nextUri = Get-SafeProperty $page '@odata.nextLink'
    }

    # Filter to members with a primary email (mail attribute) only
    $valid = [System.Collections.Generic.List[object]]::new()
    foreach ($m in $members) {
        if (Get-SafeProperty $m 'mail') {
            $valid.Add(@{
                displayName       = Get-SafeProperty $m 'displayName' '<unknown>'
                email             = (Get-SafeProperty $m 'mail').ToLower()
                userPrincipalName = Get-SafeProperty $m 'userPrincipalName'
                id                = Get-SafeProperty $m 'id'
            })
        }
        else {
            $display = Get-SafeProperty $m 'displayName' '<unknown>'
            $uid     = Get-SafeProperty $m 'id' '<unknown>'
            Write-Log WARNING "Skipping member '$display' (id: $uid) -- no primary email (mail) attribute set in Entra."
            $script:Summary.Skipped.Add("$display ($uid)")
        }
    }

    $skippedCount = $members.Count - $valid.Count
    if ($skippedCount -gt 0) {
        Write-Log WARNING "$skippedCount member(s) skipped due to missing mail attribute. They will not be synced."
    }
    if ($valid.Count -eq 0) {
        throw "No members with a primary email address found in Entra group '$EntraGroup'. Sync aborted."
    }

    Write-Log INFO "Retrieved $($valid.Count) eligible member(s) from Entra group ($($skippedCount) skipped)."
    return $valid
}

# ============================================================
# ANSYS B2C AUTHENTICATION
# ============================================================

function Get-JwtEmailClaim {
    param([string] $Token)
    try {
        $parts = $Token.Split('.')
        if ($parts.Count -lt 2) { return $null }
        $seg = $parts[1]
        $mod = $seg.Length % 4
        if ($mod -ne 0) { $seg += '=' * (4 - $mod) }
        $seg   = $seg.Replace('-', '+').Replace('_', '/')
        $bytes = [Convert]::FromBase64String($seg)
        $json  = [System.Text.Encoding]::UTF8.GetString($bytes)
        $obj   = $json | ConvertFrom-Json
        $emailProp = Get-SafeProperty $obj 'email'
        if ($emailProp) { return $emailProp }
        $upnProp = Get-SafeProperty $obj 'preferred_username'
        if ($upnProp) { return $upnProp }
        return $null
    }
    catch { return $null }
}

function Get-AuthorizationCode {
    param([string] $AuthUrl)

    # Use HttpWebRequest directly for reliable redirect and body capture.
    $location   = $null
    $body       = $null
    $statusCode = $null

    try {
        $req = [System.Net.HttpWebRequest]::Create($AuthUrl)
        $req.AllowAutoRedirect       = $false
        $req.Method                  = 'GET'
        $req.Timeout                 = 15000
        $req.AutomaticDecompression  = [System.Net.DecompressionMethods]::GZip -bor
                                       [System.Net.DecompressionMethods]::Deflate

        $resp = $null
        try {
            $resp = $req.GetResponse()
        }
        catch [System.Net.WebException] {
            $resp = $_.Exception.Response
        }

        if ($resp) {
            $statusCode = [int]$resp.StatusCode
            $location   = $resp.Headers['Location']

            # Always read the body -- needed for B2C error detail extraction
            try {
                $stream = $resp.GetResponseStream()
                if ($stream -and $stream.CanRead) {
                    $reader = [System.IO.StreamReader]::new(
                        $stream, [System.Text.Encoding]::UTF8
                    )
                    $body = $reader.ReadToEnd()
                    $reader.Close()
                }
                $stream.Close()
            }
            catch {
                Write-Log DEBUG "Could not read B2C response body: $_"
            }
            finally {
                $resp.Close()
            }
        }
    }
    catch {
        throw "Failed to reach B2C authorization endpoint: $_"
    }

    if (-not $location) {
        $detail = 'unknown'
        if ($body) {
            $match = [regex]::Match($body, '"Detail"\s*:\s*"([^"]+)"')
            if ($match.Success) { $detail = $match.Groups[1].Value }
            Write-Log ERROR "B2C returned HTTP $statusCode without redirect. Detail: $detail"
            Write-Log ERROR "B2C body excerpt: $($body.Substring(0, [Math]::Min(800, $body.Length)))"
        }
        else {
            Write-Log ERROR "B2C returned HTTP $statusCode without redirect and no readable response body."
            Write-Log ERROR "Auth URL (check for issues): $AuthUrl"
        }
        throw "B2C rejected the authorization request: $detail"
    }

    Add-Type -AssemblyName System.Web
    $uri   = [System.Uri]$location
    $query = [System.Web.HttpUtility]::ParseQueryString($uri.Query)
    $code  = $query['code']
    $err   = $query['error']

    if ($code) { return $code }
    if ($err)  {
        $desc = $query['error_description']
        throw "B2C error '$err': $desc"
    }
    throw "Redirect received but no code or error found. Location: $location"
}

function Save-TokenCache {
    param([object] $TokenResponse)
    try {
        $now    = [DateTimeOffset]::UtcNow
        # oid may not be directly in the token response; extract from id_token JWT
        $oid = Get-SafeProperty $TokenResponse 'oid' $null
        if (-not $oid -and $TokenResponse.id_token) {
            $oid = Get-JwtEmailClaim -Token $TokenResponse.id_token
            # Get-JwtEmailClaim returns email; we need oid - re-decode for it
            try {
                $parts = $TokenResponse.id_token.Split('.')
                $seg   = $parts[1]
                $mod   = $seg.Length % 4
                if ($mod -ne 0) { $seg += '=' * (4 - $mod) }
                $seg   = $seg.Replace('-', '+').Replace('_', '/')
                $obj   = [System.Text.Encoding]::UTF8.GetString([Convert]::FromBase64String($seg)) | ConvertFrom-Json
                $oidProp = $obj.PSObject.Properties['oid']
                $oid   = if ($oidProp) { $oidProp.Value } else { 'unknown' }
            }
            catch { $oid = 'unknown' }
        }
        if (-not $oid) { $oid = 'unknown' }
        $homeId = "$oid-$Policy"
        $cache  = @{ AccessToken = @{}; RefreshToken = @{}; IdToken = @{}; Account = @{} }

        if ($TokenResponse.access_token) {
            $expOn = $now.ToUnixTimeSeconds() + [int]$TokenResponse.expires_in
            $key   = "$homeId-$B2cEnv-accesstoken--$ResolvedScope--"
            $cache.AccessToken[$key] = @{
                home_account_id     = $homeId
                environment         = "$B2cEnv.b2clogin.com"
                target              = $ResolvedScope
                realm               = $TenantName
                credential_type     = 'AccessToken'
                secret              = $TokenResponse.access_token
                cached_at           = $now.ToUnixTimeSeconds().ToString()
                expires_on          = $expOn.ToString()
                extended_expires_on = ($expOn + 3600).ToString()
            }
        }
        if ($TokenResponse.refresh_token) {
            $key = "$homeId-$B2cEnv-refreshtoken-7ef9d43e-407a-40a2-94d7-14eb40416af8--"
            $cache.RefreshToken[$key] = @{
                home_account_id = $homeId
                environment     = "$B2cEnv.b2clogin.com"
                client_id       = '7ef9d43e-407a-40a2-94d7-14eb40416af8'
                credential_type = 'RefreshToken'
                secret          = $TokenResponse.refresh_token
            }
        }

        $cache | ConvertTo-Json -Depth 6 | Set-Content -Path $TokenCacheFile -Encoding UTF8
    }
    catch {
        Write-Log WARNING "Could not write token cache: $_"
    }
}

function Get-AnsysToken {
    param([string] $Pat)

    # Try token cache first
    if (Test-Path $TokenCacheFile) {
        try {
            $cache = Get-Content $TokenCacheFile -Raw | ConvertFrom-Json

            # Try access token
            $atSection = $cache.AccessToken
            if ($atSection) {
                $entry = $atSection.PSObject.Properties.Value | Select-Object -First 1
                if ($entry) {
                    $expiresOn = Get-SafeProperty $entry 'expires_on'
                    $secret = Get-SafeProperty $entry 'secret'
                    $expiry = if ($expiresOn) { [DateTimeOffset]::FromUnixTimeSeconds([long]$expiresOn) } else { [DateTimeOffset]::MinValue }
                    if ($expiry -gt [DateTimeOffset]::UtcNow.AddMinutes(5) -and $secret) {
                        $token = $secret
                        $email = Get-JwtEmailClaim -Token $token
                        Write-Log INFO "Ansys token retrieved from cache."
                        return @{ Token = $token; Email = $email }
                    }
                }
            }

            # Try refresh token
            $rtSection = $cache.RefreshToken
            if ($rtSection) {
                $rtEntry = $rtSection.PSObject.Properties.Value | Select-Object -First 1
                if ($rtEntry) {
                    Write-Log INFO "Refreshing Ansys token using cached refresh token..."
                    $tokenUri  = "https://$B2cEnv.b2clogin.com/$TenantName/$Policy/oauth2/v2.0/token"
                    $tokenBody = @{
                        grant_type    = 'refresh_token'
                        client_id     = '7ef9d43e-407a-40a2-94d7-14eb40416af8'
                        refresh_token = $rtEntry.secret
                        scope         = "openid offline_access $ResolvedScope"
                        redirect_uri  = $RedirectUri
                    }
                    try {
                        $tokenResp = Invoke-RestMethod -Uri $tokenUri -Method POST -Body $tokenBody
                        Save-TokenCache -TokenResponse $tokenResp
                        $token = Get-SafeProperty $tokenResp 'access_token'
                        $idTok = Get-SafeProperty $tokenResp 'id_token'
                        if (-not $token) { $token = $idTok }
                        if (-not $idTok) { $idTok = $token }
                        $email = Get-JwtEmailClaim -Token $idTok
                        Write-Log INFO "Ansys token refreshed and cached."
                        return @{ Token = $token; Email = $email }
                    }
                    catch {
                        Write-Log WARNING "Token refresh failed ($_) -- falling back to PAT authentication."
                    }
                }
            }
        }
        catch {
            Write-Log WARNING "Could not read token cache ($TokenCacheFile): $_"
        }
    }

    # PAT authentication
    Write-Log WARNING "No valid Ansys token cached. PAT authentication required."
    if (-not $Pat) {
        Write-Log WARNING "The Ansys PAT is a JWT token obtained from the Ansys ID Portal -- it is NOT your Ansys account password."
        $Pat = Resolve-Credential -Value '' -Prompt 'Enter your Ansys PAT' -Label 'Ansys PAT' -Secret
    }
    else {
        Write-Log INFO "Ansys PAT supplied via parameter."
    }

    $scope      = "openid offline_access $ResolvedScope"
    $scopeEnc   = [Uri]::EscapeDataString($scope)
    $redirEnc   = [Uri]::EscapeDataString($RedirectUri)
    $authUrl    = "https://$B2cEnv.b2clogin.com/$TenantName/oauth2/v2.0/authorize" +
                  "?p=$Policy" +
                  "&client_id=7ef9d43e-407a-40a2-94d7-14eb40416af8" +
                  "&nonce=defaultNonce" +
                  "&redirect_uri=$redirEnc" +
                  "&scope=$scopeEnc" +
                  "&response_type=code" +
                  "&id_token_hint=$Pat"

    Write-Log INFO "Obtaining B2C authorization code..."
    $code = Get-AuthorizationCode -AuthUrl $authUrl
    Write-Log INFO "Authorization code retrieved."

    $tokenUri  = "https://$B2cEnv.b2clogin.com/$TenantName/$Policy/oauth2/v2.0/token"
    $tokenBody = @{
        grant_type   = 'authorization_code'
        client_id    = '7ef9d43e-407a-40a2-94d7-14eb40416af8'
        code         = $code
        redirect_uri = $RedirectUri
        scope        = $scope
    }
    $tokenResp = Invoke-RestMethod -Uri $tokenUri -Method POST -Body $tokenBody

    $token = Get-SafeProperty $tokenResp 'access_token'
    $idTok = Get-SafeProperty $tokenResp 'id_token'
    if (-not $token) { $token = $idTok }
    if (-not $idTok) { $idTok = $token }
    if (-not $token) {
        $detail = $tokenResp | ConvertTo-Json -Compress
        throw "Ansys token exchange failed: $detail"
    }

    Save-TokenCache -TokenResponse $tokenResp
    Write-Log INFO "Ansys access token acquired and cached."
    $email = Get-JwtEmailClaim -Token $idTok
    return @{ Token = $token; Email = $email }
}

# ============================================================
# ANSYS IAM API HELPERS
# ============================================================

function Invoke-AnsysGet {
    param([string] $Path, [hashtable] $Query = @{})
    $uri = $AnsysIamBase + $Path
    if ($Query.Count -gt 0) {
        $qs  = ($Query.GetEnumerator() | ForEach-Object {
            '{0}={1}' -f $_.Key, [Uri]::EscapeDataString([string]$_.Value)
        }) -join '&'
        $uri = $uri + "?" + $qs
    }
    return Invoke-RestMethod -Uri $uri -Headers $script:AnsysHeaders -Method GET
}

function Invoke-AnsysPost {
    param([string] $Path, [hashtable] $Query = @{}, [object] $Body = $null)
    $uri = $AnsysIamBase + $Path
    if ($Query.Count -gt 0) {
        $qs  = ($Query.GetEnumerator() | ForEach-Object {
            '{0}={1}' -f $_.Key, [Uri]::EscapeDataString([string]$_.Value)
        }) -join '&'
        $uri = $uri + "?" + $qs
    }
    if ($null -eq $Body) {
        return Invoke-RestMethod -Uri $uri -Headers $script:AnsysHeaders -Method POST
    }
    return Invoke-RestMethod -Uri $uri -Headers $script:AnsysHeaders -Method POST `
        -Body (ConvertTo-Json -InputObject $Body -Depth 5 -Compress) -ContentType 'application/json'
}

function Invoke-AnsysPut {
    param([string] $Path, [hashtable] $Query = @{}, [object] $Body = $null)
    $uri = $AnsysIamBase + $Path
    if ($Query.Count -gt 0) {
        $qs  = ($Query.GetEnumerator() | ForEach-Object {
            '{0}={1}' -f $_.Key, [Uri]::EscapeDataString([string]$_.Value)
        }) -join '&'
        $uri = $uri + "?" + $qs
    }
    if ($DryRun) {
        $bodyStr = if ($Body) { $Body | ConvertTo-Json -Compress } else { 'null' }
        Write-Log INFO "[DRY-RUN] PUT $uri  body=$bodyStr"
        return @{}
    }
    return Invoke-RestMethod -Uri $uri -Headers $script:AnsysHeaders -Method PUT `
        -Body (ConvertTo-Json -InputObject $Body -Depth 5 -Compress) -ContentType 'application/json'
}

function Invoke-AnsysDelete {
    param([string] $Path, [object] $Body = $null)
    $uri = $AnsysIamBase + $Path
    if ($DryRun) {
        $bodyStr = if ($Body) { $Body | ConvertTo-Json -Compress } else { 'null' }
        Write-Log INFO "[DRY-RUN] DELETE $uri  body=$bodyStr"
        return
    }
    if ($Body) {
        Invoke-RestMethod -Uri $uri -Headers $script:AnsysHeaders -Method DELETE `
            -Body (ConvertTo-Json -InputObject $Body -Depth 5 -Compress) -ContentType 'application/json' | Out-Null
    }
    else {
        Invoke-RestMethod -Uri $uri -Headers $script:AnsysHeaders -Method DELETE | Out-Null
    }
}

# ============================================================
# ANSYS ID PORTAL -- ACCOUNT OPERATIONS
# ============================================================

function Find-AccountUuid {
    param([string] $CallerEmail)

    if (-not $CallerEmail) {
        throw "Cannot resolve account UUID -- caller email could not be extracted from the token. Ensure the token contains an 'email' or 'preferred_username' claim."
    }

    Write-Log INFO "Resolving account number '$AccountNumber' via user details for $CallerEmail"
    $data     = Invoke-AnsysGet -Path '/User/Details' -Query @{ email = $CallerEmail }
    $accounts = @(Get-SafeProperty $data 'accounts' @())

    foreach ($acct in $accounts) {
        $acctNum  = Get-SafeProperty $acct 'accountNumber'
        $acctName = Get-SafeProperty $acct 'accountName'
        $acctId   = Get-SafeProperty $acct 'accountId'
        if ($acctNum -eq $AccountNumber -or $acctName -eq $AccountNumber) {
            $uid = $acctId
            Write-Log INFO "Resolved account '$AccountNumber' -> UUID $uid"
            return $uid
        }
    }

    throw "Account '$AccountNumber' not found in the account list for $CallerEmail. Ensure the account exists and the PAT owner is a member of it."
}

function Get-AccountUsers {
    param([string] $AccountUuid)

    $allUsers = [System.Collections.Generic.List[object]]::new()
    $page     = 0

    do {
        $data  = Invoke-AnsysPost -Path "/Account/$AccountUuid/Users" `
                     -Query @{ PageNumber = $page; PageSize = 200 }
        $items = @(Get-SafeProperty $data 'users' @())
        foreach ($u in $items) { $allUsers.Add($u) }
        $tcPropU = Get-SafeProperty $data 'totalCount'
        $total   = if ($null -ne $tcPropU) { [int]$tcPropU } else { $allUsers.Count }
        $page++
    } while ($allUsers.Count -lt $total -and $items.Count -gt 0)

    Write-Log INFO "Account has $($allUsers.Count) current user(s)."
    return $allUsers
}

function Add-AccountUsers {
    param([string] $AccountUuid, [string[]] $Emails)

    if (-not $Emails -or $Emails.Count -eq 0) { return }
    Write-Log INFO "Bulk importing $($Emails.Count) user(s) into account $AccountUuid"

    $payload = @($Emails | ForEach-Object { @{ email = $_ } })
    $qs      = 'removeOthers=false'
    if (-not $SendEmail) { $qs += '&SuppressEmail=true' }

    if ($DryRun) {
        Write-Log INFO "[DRY-RUN] PUT /Account/$AccountUuid/Users?$qs  first 3: $(($payload | Select-Object -First 3 | ConvertTo-Json -Compress))"
        return
    }

    $uri    = "$AnsysIamBase/Account/$AccountUuid/Users?$qs"
    $result = Invoke-RestMethod -Uri $uri -Headers $script:AnsysHeaders -Method PUT `
                  -Body (ConvertTo-Json -InputObject $payload -Depth 3 -Compress) -ContentType 'application/json'

    foreach ($bucket in @('created','updated','rejected','failed')) {
        $entries = Get-SafeProperty $result $bucket @()
        if ($entries -and @($entries).Count -gt 0) {
            Write-Log INFO "  Account import $bucket`: $(@($entries).Count) user(s)"
            foreach ($e in $entries) {
                $errMsg = Get-SafeProperty $e 'errorMessage'
                if ($errMsg) { Write-Log WARNING "    $(Get-SafeProperty $e 'email') -- $errMsg" }
            }
        }
    }
}

function Remove-AccountUsers {
    param([string] $AccountUuid, [string[]] $UserUuids)

    if (-not $UserUuids -or $UserUuids.Count -eq 0) { return }
    Write-Log INFO "  - Removing $($UserUuids.Count) user(s) from account $AccountUuid"
    Invoke-AnsysDelete -Path "/Account/$AccountUuid/Users" -Body @($UserUuids)
}

# ============================================================
# ANSYS ID PORTAL -- GROUP OPERATIONS
# ============================================================

function Find-AnsysGroup {
    param([string] $AccountUuid, [string] $GrpName)

    $data  = Invoke-AnsysPost -Path "/Account/$AccountUuid/Groups" `
                 -Query @{ searchValue = $GrpName; PageSize = 50; PageNumber = 0 }
    $items = @(Get-SafeProperty $data 'groups' @())

    foreach ($g in $items) {
        if ((Get-SafeProperty $g 'name') -eq $GrpName) {
            Write-Log INFO "Found group '$GrpName' -> UUID $(Get-SafeProperty $g 'id')"
            return $g
        }
    }
    return $null
}

function New-AnsysGroup {
    param([string] $AccountUuid, [string] $GrpName)

    Write-Log INFO "Creating group '$GrpName' under account $AccountUuid"
    $body = @{
        name        = $GrpName
        description = "Synced from Entra ID group: $GrpName"
        groupType   = 1
    }
    if ($DryRun) {
        Write-Log INFO "[DRY-RUN] POST /Account/$AccountUuid/Group  body=$($body | ConvertTo-Json -Compress)"
        return @{ groupId = 'dry-run-group-uuid' }
    }
    $result = Invoke-AnsysPost -Path "/Account/$AccountUuid/Group" -Body $body
    Write-Log INFO "Created group '$GrpName' -> UUID $(Get-SafeProperty $result 'groupId')"
    return $result
}

function Get-GroupMembers {
    param([string] $AccountUuid, [string] $GroupUuid)

    $allMembers = [System.Collections.Generic.List[object]]::new()
    $page       = 0

    do {
        $data  = Invoke-AnsysPost -Path "/Account/$AccountUuid/Group/$GroupUuid/Members" `
                     -Query @{ PageNumber = $page; PageSize = 200 }
        $items = @(Get-SafeProperty $data 'members' @())
        foreach ($m in $items) { $allMembers.Add($m) }
        $tcPropM = Get-SafeProperty $data 'totalCount'
        $total   = if ($null -ne $tcPropM) { [int]$tcPropM } else { $allMembers.Count }
        $page++
    } while ($allMembers.Count -lt $total -and $items.Count -gt 0)

    Write-Log INFO "Group has $($allMembers.Count) current member(s)."
    return $allMembers
}

function Add-GroupMembers {
    param([string] $AccountUuid, [string] $GroupUuid, [string[]] $Emails)

    if (-not $Emails -or $Emails.Count -eq 0) { return }
    Write-Log INFO "Importing $($Emails.Count) user(s) into group $GroupUuid"

    $body = @{ users = @($Emails | ForEach-Object { @{ email = $_; isAdmin = $false } }) }

    if ($DryRun) {
        Write-Log INFO "[DRY-RUN] PUT /Account/$AccountUuid/Group/$GroupUuid/Members/Import  body=$($body | ConvertTo-Json -Compress)"
        return
    }

    $uri    = "$AnsysIamBase/Account/$AccountUuid/Group/$GroupUuid/Members/Import?removeOthers=false"
    $result = Invoke-RestMethod -Uri $uri -Headers $script:AnsysHeaders -Method PUT `
                  -Body (ConvertTo-Json -InputObject $body -Depth 4 -Compress) -ContentType 'application/json'

    foreach ($bucket in @('success','notModified','error')) {
        $entries = Get-SafeProperty $result $bucket @()
        if ($entries -and @($entries).Count -gt 0) {
            Write-Log INFO "  Group import $bucket`: $(@($entries).Count) user(s)"
            foreach ($e in $entries) {
                $errMsg = Get-SafeProperty $e 'errorMessage'
                if ($errMsg) { Write-Log WARNING "    $(Get-SafeProperty $e 'email') -- $errMsg" }
            }
        }
    }
}

function Remove-GroupMembers {
    param([string] $AccountUuid, [string] $GroupUuid, [string[]] $MemberUuids)

    if (-not $MemberUuids -or $MemberUuids.Count -eq 0) { return }
    Write-Log INFO "  - Removing $($MemberUuids.Count) member(s) from group $GroupUuid"
    $body = @{ members = @($MemberUuids | ForEach-Object { @{ memberId = $_ } }) }
    Invoke-AnsysDelete -Path "/Account/$AccountUuid/Group/$GroupUuid/Members" -Body $body
}

# ============================================================
# SYNC -- ACCOUNT
# ============================================================

function Sync-ToAccount {
    param([string] $AccountUuid, [object[]] $EntraMembers)

    Write-Log INFO "=== Syncing to ACCOUNT $AccountUuid ==="

    $desiredEmails = [System.Collections.Generic.HashSet[string]]::new(
        [System.StringComparer]::OrdinalIgnoreCase
    )
    foreach ($m in $EntraMembers) { $desiredEmails.Add($m.email) | Out-Null }
    Write-Log INFO "Desired: $($desiredEmails.Count) user(s)"

    $currentMembers = Get-AccountUsers -AccountUuid $AccountUuid
    $currentByEmail = @{}
    foreach ($m in $currentMembers) {
        $mEmail  = Get-SafeProperty $m 'email'
        $mUserId = Get-SafeProperty $m 'userId'
        if ($mEmail -and $mUserId) {
            $currentByEmail[$mEmail.ToLower()] = $mUserId
        }
    }

    $toAdd    = @($desiredEmails | Where-Object { -not $currentByEmail.ContainsKey($_) })
    $toRemove = @($currentByEmail.Keys | Where-Object { -not $desiredEmails.Contains($_) })

    if ($script:PreserveEmails.Count -gt 0) {
        $protected = @($toRemove | Where-Object { $script:PreserveEmails.Contains($_) })
        if ($protected.Count -gt 0) {
            Write-Log INFO "Skipping removal of $($protected.Count) preserved address(es): $($protected -join ', ')"
        }
        $toRemove = @($toRemove | Where-Object { -not $script:PreserveEmails.Contains($_) })
    }

    Write-Log INFO "To add: $($toAdd.Count),  To remove: $($toRemove.Count)"
    foreach ($e in ($toAdd    | Sort-Object)) { $script:Summary.Added.Add($e) }
    foreach ($e in ($toRemove | Sort-Object)) { $script:Summary.Removed.Add($e) }

    Add-AccountUsers -AccountUuid $AccountUuid -Emails $toAdd

    if ($toRemove.Count -gt 0) {
        $uuidsToRemove = @($toRemove | ForEach-Object { $currentByEmail[$_] })
        Remove-AccountUsers -AccountUuid $AccountUuid -UserUuids $uuidsToRemove
    }

    Write-Log INFO "Account sync complete."
}

# ============================================================
# SYNC -- GROUP
# ============================================================

function Sync-ToGroup {
    param([string] $AccountUuid, [string] $GrpName, [object[]] $EntraMembers)

    Write-Log INFO "=== Syncing to GROUP '$GrpName' under account $AccountUuid ==="

    $group = Find-AnsysGroup -AccountUuid $AccountUuid -GrpName $GrpName
    if (-not $group) {
        Write-Log INFO "Group '$GrpName' not found -- creating it."
        $group = New-AnsysGroup -AccountUuid $AccountUuid -GrpName $GrpName
    }

    # find returns id; create returns groupId
    $groupUuid = Get-SafeProperty $group 'groupId'
    if (-not $groupUuid) { $groupUuid = Get-SafeProperty $group 'id' }

    if (-not $groupUuid -or $groupUuid -eq 'dry-run-group-uuid') {
        if ($DryRun) {
            Write-Log INFO "[DRY-RUN] Skipping member sync -- no real group UUID."
            return
        }
        throw "Could not resolve group UUID. Response: $($group | ConvertTo-Json -Compress)"
    }
    Write-Log INFO "Using group UUID: $groupUuid"

    $desiredEmails = [System.Collections.Generic.HashSet[string]]::new(
        [System.StringComparer]::OrdinalIgnoreCase
    )
    foreach ($m in $EntraMembers) { $desiredEmails.Add($m.email) | Out-Null }
    Write-Log INFO "Desired: $($desiredEmails.Count) member(s)"

    $currentMembers = Get-GroupMembers -AccountUuid $AccountUuid -GroupUuid $groupUuid
    $currentByEmail = @{}
    foreach ($m in $currentMembers) {
        if (Get-SafeProperty $m 'isGroupAsMember') { continue }
        $gmEmail    = Get-SafeProperty $m 'email'
        $gmMemberId = Get-SafeProperty $m 'memberId'
        if ($gmEmail -and $gmMemberId) {
            $currentByEmail[$gmEmail.ToLower()] = $gmMemberId
        }
    }

    $toAdd    = @($desiredEmails | Where-Object { -not $currentByEmail.ContainsKey($_) })
    $toRemove = @($currentByEmail.Keys | Where-Object { -not $desiredEmails.Contains($_) })

    if ($script:PreserveEmails.Count -gt 0) {
        $protected = @($toRemove | Where-Object { $script:PreserveEmails.Contains($_) })
        if ($protected.Count -gt 0) {
            Write-Log INFO "Skipping removal of $($protected.Count) preserved address(es): $($protected -join ', ')"
        }
        $toRemove = @($toRemove | Where-Object { -not $script:PreserveEmails.Contains($_) })
    }

    Write-Log INFO "To add: $($toAdd.Count),  To remove: $($toRemove.Count)"
    foreach ($e in ($toAdd    | Sort-Object)) { $script:Summary.Added.Add($e) }
    foreach ($e in ($toRemove | Sort-Object)) { $script:Summary.Removed.Add($e) }

    Add-GroupMembers -AccountUuid $AccountUuid -GroupUuid $groupUuid -Emails $toAdd

    if ($toRemove.Count -gt 0) {
        $uuidsToRemove = @($toRemove | ForEach-Object { $currentByEmail[$_] })
        Remove-GroupMembers -AccountUuid $AccountUuid -GroupUuid $groupUuid -MemberUuids $uuidsToRemove
    }

    Write-Log INFO "Group sync complete."
}

# ============================================================
# MAIN
# ============================================================

# Validate arguments
if ($TargetType -eq 'group' -and -not $GroupName) {
    Write-Log ERROR "-GroupName is required when -TargetType is 'group'."
    exit 1
}
if (($LogRecipient -or $WarnRecipient) -and -not $SmtpRelay) {
    Write-Log WARNING "Email recipients specified but -SmtpRelay not provided. No emails will be sent."
}
if (($LogRecipient -or $WarnRecipient) -and -not $SmtpFrom) {
    Write-Log WARNING "Email recipients specified but -SmtpFrom not provided. No emails will be sent."
}

Add-Type -AssemblyName System.Web
Open-LogFile -Path $LogFile

$script:Summary.StartedAt = (Get-Date).ToUniversalTime().ToString('o')

$logRecipients  = Split-Recipients -Value $LogRecipient
$warnRecipients = Split-Recipients -Value $WarnRecipient

# Build preserve email set
$script:PreserveEmails = [System.Collections.Generic.HashSet[string]]::new(
    [System.StringComparer]::OrdinalIgnoreCase
)
if ($PreserveEmail) {
    foreach ($e in ($PreserveEmail -split ',' | ForEach-Object { $_.Trim() } | Where-Object { $_ })) {
        $script:PreserveEmails.Add($e.ToLower()) | Out-Null
    }
    Write-Log INFO "Preserved emails (never removed): $(($script:PreserveEmails | Sort-Object) -join ', ')"
}

try {
    # 1. Resolve Entra credentials
    $EntraClientId     = Resolve-Credential -Value $EntraClientId -Prompt 'Enter Entra App Registration Client ID' -Label 'Entra client ID'
    $EntraClientSecret = Resolve-Credential -Value $EntraClientSecret -Prompt 'Enter Entra App Registration Client Secret' -Label 'Entra client secret' -Secret

    # 2. Fetch Entra group members
    $graphToken   = Get-GraphToken -ClientId $EntraClientId -ClientSecret $EntraClientSecret
    $entraMembers = @(Get-EntraGroupMembers -GraphToken $graphToken)

    if ($entraMembers.Count -eq 0) {
        Write-Log WARNING "Entra group returned 0 eligible members. Nothing to sync."
        $script:Summary.ExitStatus = 'no-op'
        return
    }

    # 3. Obtain Ansys token
    $tokenResult = Get-AnsysToken -Pat $AnsysPat
    $ansysToken  = $tokenResult.Token
    $callerEmail = $tokenResult.Email

    if ($callerEmail) {
        Write-Log INFO "Authenticated as: $callerEmail"
    }
    else {
        Write-Log WARNING "Could not extract email from token -- account UUID resolution will fail."
    }

    # Set Ansys API headers
    $script:AnsysHeaders = @{
        Authorization = "Bearer $ansysToken"
        Accept        = 'application/json'
    }

    # 4. Resolve account UUID
    $accountUuid = Find-AccountUuid -CallerEmail $callerEmail

    # 5. Sync
    if ($TargetType -eq 'account') {
        Sync-ToAccount -AccountUuid $accountUuid -EntraMembers $entraMembers
    }
    else {
        Sync-ToGroup -AccountUuid $accountUuid -GrpName $GroupName -EntraMembers $entraMembers
    }
}
catch {
    $err = "Unhandled exception: $_"
    Write-Log ERROR $err
    $script:Summary.Errors.Add($err)
    $script:Summary.ExitStatus = 'failure'
}
finally {
    $script:Summary.FinishedAt = (Get-Date).ToUniversalTime().ToString('o')
    Write-Log INFO ('Sync finished -- status: {0}  added: {1}  removed: {2}  skipped: {3}  errors: {4}' -f
        $script:Summary.ExitStatus,
        $script:Summary.Added.Count,
        $script:Summary.Removed.Count,
        $script:Summary.Skipped.Count,
        $script:Summary.Errors.Count
    )
    Send-WarnEmail -Recipients $warnRecipients
    Send-LogEmail  -Recipients $logRecipients
    Close-LogFile
}

if ($script:Summary.ExitStatus -eq 'failure') { exit 1 }