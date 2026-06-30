# HapDup External Submodules

The original HapDup project uses Git submodules for external tools:

- Flye: https://github.com/fenderglass/Flye
- Margin: https://github.com/UCSC-nanopore-cgl/margin
- PEPPER: https://github.com/kishwarshafin/pepper

This BayesPhase repository keeps the HapDup integration code and `.gitmodules` metadata, but does not vendor the full third-party tool source trees. Install or initialize these tools from their upstream projects when running the HapDup-BayesPhase workflow.
