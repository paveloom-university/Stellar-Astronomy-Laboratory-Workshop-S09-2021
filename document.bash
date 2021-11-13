#!/bin/bash

# A script to document the crate

cargo doc
RUSTDOCFLAGS="--html-in-header assets/katex-header.html" cargo doc --no-deps
echo "<meta http-equiv=\"refresh\" content=\"0; url=kidou\">" > target/doc/index.html
