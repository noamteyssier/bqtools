precommit: check check-all clippy fmt test test-all

check:
    cargo check

check-all:
    cargo check --all-features

clippy:
    cargo clippy --all-features -- -D warnings

fmt:
    cargo fmt

test:
    cargo test

test-all:
    cargo test --all-features

install:
    cargo install --path .

docs:
    cargo run --all-features --bin bqdoc
    uvx zensical build
    rsync -a --include='*/' --include='*.md' --exclude='*' docs/ site/

docs-serve:
    cargo run --all-features --bin bqdoc
    uvx zensical serve
