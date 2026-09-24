use std::process::{Child, Command};

use anyhow::Result;
use log::warn;

use crate::cli::FileFormat;

use super::{utils::name_fifo, PairedChannels, RecordPair};

#[derive(Clone, Copy)]
pub enum ExecMode<'a> {
    /// `-x`: one shell invocation per FIFO (or per R1/R2 pair for paired files).
    PerFifo(&'a str),
    /// `-X`: one shell invocation with all FIFO paths substituted in.
    Batch(&'a str),
}

/// Validate that the exec template contains the substitution tokens required for
/// the file type, so a missing `{}` / `{R1}` / `{R2}` fails fast rather than
/// leaving a FIFO open with no reader (which hangs indefinitely).
///
/// For paired files, at least one of `{R1}` or `{R2}` is required. Supplying
/// only one is valid — only that channel's FIFO will be created and written.
pub fn validate_template(template: &str, paired: bool) -> Result<()> {
    if paired {
        if !template.contains("{R1}") && !template.contains("{R2}") {
            anyhow::bail!(
                "exec template for a paired file must contain at least one of {{R1}} or {{R2}}"
            );
        }
    } else if !template.contains("{}") {
        anyhow::bail!("exec template must contain {{}} for the FIFO path");
    }
    Ok(())
}

/// Returns which paired channels the template requires.
pub fn required_channels(template: &str) -> PairedChannels {
    match (template.contains("{R1}"), template.contains("{R2}")) {
        (true, true) => PairedChannels::Both,
        (true, false) => PairedChannels::R1Only,
        (false, true) => PairedChannels::R2Only,
        (false, false) => unreachable!(
            "exec template for a single file must contain at least one of {{R1}} or {{R2}}"
        ),
    }
}

/// Spawn consumer subprocesses according to `mode`, returning their handles.
///
/// Must be called after FIFOs are created but before writer threads are spawned,
/// because opening a FIFO for writing blocks until a reader connects.
pub fn spawn_consumers(
    mode: ExecMode<'_>,
    basename: &str,
    paired: bool,
    num_pipes: usize,
    format: FileFormat,
) -> Result<Vec<Child>> {
    let mut children = Vec::new();
    match mode {
        ExecMode::PerFifo(template) => {
            for pid in 0..num_pipes {
                let cmd = if paired {
                    let r1 = name_fifo(basename, pid, RecordPair::R1, format);
                    let r2 = name_fifo(basename, pid, RecordPair::R2, format);
                    template.replace("{R1}", &r1).replace("{R2}", &r2)
                } else {
                    let path = name_fifo(basename, pid, RecordPair::Unpaired, format);
                    template.replace("{}", &path)
                }
                .replace("{n}", &pid.to_string());
                children.push(sh(&cmd)?);
            }
        }
        ExecMode::Batch(template) => {
            if template.contains("{n}") {
                warn!(
                    "{{n}} is not expanded by --exec-batch; did you mean -x/--exec (one command per pipe)?"
                );
            }
            let cmd = if paired {
                let r1s: Vec<_> = (0..num_pipes)
                    .map(|pid| name_fifo(basename, pid, RecordPair::R1, format))
                    .collect();
                let r2s: Vec<_> = (0..num_pipes)
                    .map(|pid| name_fifo(basename, pid, RecordPair::R2, format))
                    .collect();
                // An adjacent `{R1} {R2}` expands as interleaved pairs (r1_0 r2_0 r1_1 r2_1 …)
                // so positional paired-argument tools receive each pair together; any
                // remaining `{R1}` / `{R2}` expand to their own space-joined lists.
                let interleaved: Vec<_> = r1s
                    .iter()
                    .zip(&r2s)
                    .flat_map(|(r1, r2)| [r1.as_str(), r2.as_str()])
                    .collect();
                template
                    .replace("{R1} {R2}", &interleaved.join(" "))
                    .replace("{R1}", &r1s.join(" "))
                    .replace("{R2}", &r2s.join(" "))
            } else {
                let paths: Vec<_> = (0..num_pipes)
                    .map(|pid| name_fifo(basename, pid, RecordPair::Unpaired, format))
                    .collect();
                template.replace("{}", &paths.join(" "))
            };
            children.push(sh(&cmd)?);
        }
    }
    Ok(children)
}

fn sh(cmd: &str) -> Result<Child> {
    log::debug!("exec: sh -c {cmd:?}");
    Command::new("sh")
        .arg("-c")
        .arg(cmd)
        .spawn()
        .map_err(Into::into)
}
