#!/usr/bin/env python3
import os
import time

import requests


BASE_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
PMC_IDCONV = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"
NCBI_EMAIL = os.getenv("NCBI_EMAIL", "your.name@institution.edu")
NCBI_KEY = os.getenv("NCCI_API_KEY") or os.getenv("NCBI_API_KEY")
NCBI_MAX_RETRIES = int(os.getenv("NCBI_MAX_RETRIES", "5"))
NCBI_RETRY_STATUS_CODES = {429, 500, 502, 503, 504}

# Be conservative here. These scripts make many sequential metadata + efetch calls.
DEFAULT_MIN_INTERVAL = "0.75" if NCBI_KEY else "1.0"
NCBI_MIN_INTERVAL_SECONDS = float(os.getenv("NCBI_MIN_INTERVAL_SECONDS", DEFAULT_MIN_INTERVAL))


class NCBIClient:
    def __init__(self, tool_name):
        self.tool_name = tool_name
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": f"{tool_name} ({NCBI_EMAIL})"})
        self._last_request_at = 0.0

    def default_params(self, extra=None):
        params = {"tool": self.tool_name, "email": NCBI_EMAIL}
        if NCBI_KEY:
            params["api_key"] = NCBI_KEY
        if extra:
            params.update(extra)
        return params

    def _throttle(self):
        elapsed = time.monotonic() - self._last_request_at
        remaining = NCBI_MIN_INTERVAL_SECONDS - elapsed
        if remaining > 0:
            time.sleep(remaining)

    def request(self, method, url, *, params=None, data=None, timeout=30, retries=None):
        retries = NCBI_MAX_RETRIES if retries is None else retries
        last_exc = None

        for attempt in range(1, retries + 1):
            self._throttle()
            try:
                resp = self.session.request(method, url, params=params, data=data, timeout=timeout)
                self._last_request_at = time.monotonic()
                if resp.status_code in NCBI_RETRY_STATUS_CODES:
                    raise requests.HTTPError(
                        f"{resp.status_code} Server Error for url: {resp.url}",
                        response=resp,
                    )
                resp.raise_for_status()
                return resp
            except (requests.Timeout, requests.ConnectionError, requests.HTTPError) as exc:
                self._last_request_at = time.monotonic()
                last_exc = exc
                status_code = getattr(getattr(exc, "response", None), "status_code", None)
                retryable = status_code in NCBI_RETRY_STATUS_CODES or isinstance(
                    exc,
                    (requests.Timeout, requests.ConnectionError),
                )
                if not retryable or attempt == retries:
                    raise
                wait_seconds = min(max(2 ** attempt, NCBI_MIN_INTERVAL_SECONDS), 30.0)
                print(
                    f"  [NCBI retry {attempt}/{retries}] {method.upper()} failed ({exc}). "
                    f"Sleeping {wait_seconds:.1f}s before retry."
                )
                time.sleep(wait_seconds)

        raise last_exc

    def get(self, url, *, params=None, timeout=30, retries=None):
        return self.request("get", url, params=params, timeout=timeout, retries=retries)

    def post(self, url, *, data=None, timeout=30, retries=None):
        return self.request("post", url, data=data, timeout=timeout, retries=retries)
