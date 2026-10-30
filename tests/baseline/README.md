# Frozen test baselines

These two source files are test fixtures from the pre-mixed implementations.
They are never imported by the production entry points. Exact-default parity
checks compare every preexisting result against these fixtures.

Canonical UTF-8/LF SHA-256 values (normalize CRLF to LF before hashing):

- `R/summclust.R`: `99661889d38fa692b0726d02ba9464fa0cb50834d48c98684d869ca0a9f7c41c`.
- `python/summclust.py`: `bf6f5a111a573ffbad01e9ff5fce1f184adf87d702a9d547b49088ae2e011e5f`.

The originals were standalone base-R and NumPy ports of Stata summclust 5.1.1
with the September 30, 2026 presentation revision, frozen before adding CV31.
The same repository license and source attribution apply.
