//! Generates the `MkDocs` command reference from the clap definitions.
//!
//! Run from the repo root: `cargo run --features bqdoc --bin bqdoc`
use std::fmt::Write as _;
use std::fs;
use std::path::Path;

use anyhow::Result;
use bqtools::cli::Cli;
use clap::{Arg, ArgAction, CommandFactory};

fn takes_value(arg: &Arg) -> bool {
    !matches!(
        arg.get_action(),
        ArgAction::SetTrue | ArgAction::SetFalse | ArgAction::Count | ArgAction::Help
    )
}

fn esc(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

/// Escaped text with `code` spans and blank-line separated paragraphs as HTML.
fn html(text: &str) -> String {
    let mut out = String::new();
    for para in text.split("\n\n") {
        let para = para.split_whitespace().collect::<Vec<_>>().join(" ");
        out.push_str("<p>");
        for (i, t) in esc(&para).split('`').enumerate() {
            if i % 2 == 1 {
                let _ = write!(out, "<code>{t}</code>");
            } else {
                out.push_str(t);
            }
        }
        out.push_str("</p>");
    }
    out
}

fn term(arg: &Arg) -> String {
    let mut names = Vec::new();
    if let Some(l) = arg.get_long() {
        names.push(format!("<code>--{l}</code>"));
    }
    if let Some(s) = arg.get_short() {
        names.push(format!("<code>-{s}</code>"));
    }
    let value = if takes_value(arg) {
        arg.get_value_names().map_or_else(
            || arg.get_id().as_str().to_uppercase(),
            |v| {
                v.iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>()
                    .join(" ")
            },
        )
    } else {
        String::new()
    };
    if names.is_empty() {
        return format!("<code>{}</code>", esc(&value));
    }
    let names = names.join(", ");
    if value.is_empty() {
        names
    } else {
        format!(
            "{names} <i>{}</i>",
            esc(&value.to_lowercase().replace('_', "-"))
        )
    }
}

fn write_args(out: &mut String, cmd: &str, args: &[&Arg]) {
    out.push_str("<dl class=\"cli-reference\">\n");
    for arg in args {
        let id = format!("{cmd}--{}", arg.get_id().as_str().replace('_', "-"));
        let _ = write!(
            out,
            "<dt id=\"{id}\"><a href=\"#{id}\">{}</a></dt>\n<dd>",
            term(arg)
        );
        if let Some(help) = arg.get_long_help().or_else(|| arg.get_help()) {
            out.push_str(&html(&help.to_string()));
        }
        let defaults: Vec<_> = arg
            .get_default_values()
            .iter()
            .map(|d| d.to_string_lossy().into_owned())
            .collect();
        if takes_value(arg) && !defaults.is_empty() {
            let _ = write!(
                out,
                "<p>Default: <code>{}</code></p>",
                esc(&defaults.join(","))
            );
        }
        let possible: Vec<_> = arg
            .get_possible_values()
            .into_iter()
            .filter(|p| !p.is_hide_set())
            .collect();
        if !possible.is_empty() {
            out.push_str("<p>Possible values:</p><ul>");
            for p in possible {
                let _ = write!(out, "<li><code>{}</code>", esc(p.get_name()));
                if let Some(h) = p.get_help() {
                    let _ = write!(out, ": {}", esc(&h.to_string()));
                }
                out.push_str("</li>");
            }
            out.push_str("</ul>");
        }
        out.push_str("</dd>\n");
    }
    out.push_str("</dl>\n");
}

fn page(cmd: &mut clap::Command) -> String {
    let name = cmd.get_name().to_string();
    let mut out = format!("# bqtools {name}\n\n");
    if let Some(about) = cmd.get_long_about().or_else(|| cmd.get_about()) {
        let _ = write!(out, "{about}\n\n");
    }
    let usage = cmd.render_usage().to_string();
    let _ = write!(
        out,
        "## Usage\n\n```text\n{}\n```\n",
        usage.trim_start_matches("Usage: ")
    );

    let args: Vec<&Arg> = cmd
        .get_arguments()
        .filter(|a| !matches!(a.get_id().as_str(), "help" | "version") && !a.is_hide_set())
        .collect();
    let mut headings: Vec<Option<&str>> = Vec::new();
    for a in &args {
        let h = a.get_help_heading();
        if !headings.contains(&h) {
            headings.push(h);
        }
    }
    for h in headings {
        let group: Vec<&Arg> = args
            .iter()
            .copied()
            .filter(|a| a.get_help_heading() == h)
            .collect();
        let title = h.map_or_else(
            || "Arguments".to_string(),
            |h| {
                let mut c = h.to_lowercase();
                c[..1].make_ascii_uppercase();
                c
            },
        );
        let _ = write!(out, "\n## {title}\n\n");
        write_args(&mut out, &name, &group);
    }
    out
}

fn main() -> Result<()> {
    let mut cli = Cli::command();
    cli.build();

    let dir = Path::new("docs/commands");
    fs::create_dir_all(dir)?;

    // The nav lives in mkdocs.yml; fail if a command is missing from it.
    let nav = fs::read_to_string("mkdocs.yml")?;
    let names: Vec<String> = cli
        .get_subcommands()
        .filter(|c| !c.is_hide_set() && c.get_name() != "help")
        .map(|c| c.get_name().to_string())
        .collect();
    for name in names {
        let sub = cli.find_subcommand_mut(&name).expect("subcommand exists");
        fs::write(dir.join(format!("{name}.md")), page(sub))?;
        anyhow::ensure!(
            nav.contains(&format!("commands/{name}.md")),
            "add commands/{name}.md to the nav in mkdocs.yml"
        );
    }
    Ok(())
}
