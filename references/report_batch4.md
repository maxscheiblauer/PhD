# Download Report — Batch 4 (VPN session, 2026-10-02)

Context: user connected via TU Wien VPN, hoping institutional IP-based access would
unlock the items earlier batches skipped as paywalled. Confirmed the VPN is active
and routing through TU Wien's network (public IP 128.130.8.38, AS679 Technische
Universitat Wien, Vienna).

## Newly downloaded (3)

- **oguamalam2024a** — was never attempted in batches 1–3 (excluded by mistake as
  "out of range"). Freely available on arXiv: https://arxiv.org/pdf/2412.19004
- **mahalanobis1936** — P. Mahalanobis, "On the generalised distance in statistics,"
  Proc. National Institute of Sciences of India, Vol. 2 (1936). Found on INSA's own
  open-access archive: https://insa.nic.in/writereaddata/UpLoadedFiles/PINSA/Vol02_1936_1_Art05.pdf
- (UVaDOC and HAL retries below did not succeed, but are included for completeness.)

## Key finding: VPN does not unlock the Cloudflare-gated publishers

Elsevier (ScienceDirect), Wiley (Online Library), Taylor & Francis (tandfonline), and
Oxford University Press (academic.oup.com) all sit behind Cloudflare's managed JS
challenge ("Just a moment..."), which returns HTTP 403 to any non-browser HTTP client
— regardless of the client's IP address or institutional network. This blocks `curl`/
scripted access even while genuinely on the TU Wien network; it would require a real,
JS-executing browser session (ideally logged into TU Wien's library proxy/Shibboleth)
to get past the challenge and then check for entitlement.

Springer (link.springer.com) is not Cloudflare-walled, so pages load via `curl`, but
Springer's access entitlement is granted through a Shibboleth/IdP login flow
(`idp.springer.com`), not plain IP recognition — confirmed by the response's
`sim-inst-token=""` (empty institutional token). So even though the Springer pages
load, no access was actually granted for `gupta2012`, `ramsay2005`, `rousseeuw1985`,
or `young1985`.

Practical implication: **no browser-automation tool is connected to this session**
(no Chrome extension / computer-use), so none of the Cloudflare- or Shibboleth-gated
items can be fetched programmatically from here, VPN or not. To get these, the fastest
path is for the user to open the TU Wien library's resource/EZproxy link (or the DOI
directly while the VPN is active) in their own browser and save the PDF manually — a
real logged-in browser session clears both the bot-challenge and the entitlement check
at once.

## Still missing (17) — attempted again this session, outcome unchanged

Paywalled, blocked by Cloudflare bot-challenge (Elsevier/Wiley/T&F/OUP) even on VPN —
**recommend manual browser download via TU Wien library access**:
- ash2014 — Topics in Stochastic Processes (Elsevier book), DOI 10.1016/C2013-0-06165-8
- dawid1981 — Biometrika/OUP, DOI 10.2307/2335827
- dutilleul1999 — J. Stat. Comput. Simul. (T&F), DOI 10.1080/00949659908811970
- jacques2014b — Comput. Stat. & Data Analysis (Elsevier), DOI 10.1016/j.csda.2012.12.004
- manceur2013 — J. Comput. Appl. Math. (Elsevier), DOI 10.1016/j.cam.2012.09.017
- maronna2019 — Robust Statistics (Wiley book), DOI 10.1002/9781119214656
- rousseeuw1999 — Technometrics (T&F), DOI 10.1080/00401706.1999.10485670
- rousseeuw2005 — Robust statistics... influence functions (Wiley book), DOI 10.1002/9781118186435
- wang2024 — WIREs Comput. Stat. (Wiley), DOI 10.1002/wics.1638
- yu2025 — WIREs Comput. Stat. (Wiley), DOI 10.1002/wics.70013

Paywalled, Springer-hosted — loads via curl but no entitlement without Shibboleth
login (same recommendation: manual browser download while on VPN, logging in via
TU Wien's SSO if prompted):
- gupta2012 — Elliptically Contoured Models in Statistics, DOI 10.1007/978-94-011-1646-6
- ramsay2005 — Functional Data Analysis, DOI 10.1007/b98888
- rousseeuw1985 — Multivariate estimation with high breakdown point, DOI 10.1007/978-94-009-5438-0_20
- young1985 — Monotonic solutions of cooperative games, DOI 10.1007/BF01769885

Open-access copies exist but the host itself is unreachable from this network (not a
paywall issue — retried this session, still fails):
- garciaescudero2024 — UVaDOC (uva.es) connection times out (TCP-level, not HTTP)
- garciaescudero2022 — same UVaDOC host
- jacques2014a — HAL (hal.science) now returns its own bot-check page ("Making sure
  you're not a bot!") rather than the Anubis block seen before — still not passable
  without a real browser

## Summary

46/63 references now saved locally (up from 44/63 before this session). The remaining
17 all require either a real browser session (for the Cloudflare/Shibboleth-gated
publishers) or a working connection to UVaDOC/HAL (transient host issues) — none of
which an automated VPN+curl session can resolve on its own.
