#!/bin/bash

# A script to document the crate

cargo doc
RUSTDOCFLAGS="--html-in-header assets/katex-header.html" cargo doc --no-deps
mv target/doc public
echo "<meta http-equiv=\"refresh\" content=\"0; url=kidou\">" >public/index.html
