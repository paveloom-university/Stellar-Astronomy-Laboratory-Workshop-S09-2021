# Notices

To build this crate's documentation with [KaTeX](https://katex.org/) support, but all other dependencies without it, run:

```bash
cargo doc
RUSTDOCFLAGS="--html-in-header assets/katex-header.html" cargo doc --no-deps --open
```
