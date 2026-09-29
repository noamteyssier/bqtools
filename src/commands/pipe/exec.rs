use std::process::{Child, Command};

use anyhow::Result;
use log::warn;

use crate::cli::FileFormat;

use super::{utils::name_fifo, PairedChannels, RecordPair};

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

/// Spawn consumer subprocesses, returning their handles.
///
/// `batch` selects `-X` (one shell invocation with all FIFO paths substituted in)
/// over `-x` (one invocation per FIFO, or per R1/R2 pair for paired files).
///
/// Must be called after FIFOs are created but before writer threads are spawned,
/// because opening a FIFO for writing blocks until a reader connects.
pub fn spawn_consumers(
    template: &str,
    batch: bool,
    basename: &str,
    paired: bool,
    num_pipes: usize,
    format: FileFormat,
) -> Result<Vec<Child>> {
    let paths = |pair| -> Vec<String> {
        (0..num_pipes)
            .map(|pid| name_fifo(basename, pid, pair, format))
            .collect()
    };
    let sh = |cmd: &str| -> Result<Child> {
        log::debug!("exec: sh -c {cmd:?}");
        Ok(Command::new("sh").arg("-c").arg(cmd).spawn()?)
    };
    if !batch {
        let (a, b) = if paired {
            (paths(RecordPair::R1), paths(RecordPair::R2))
        } else {
            (paths(RecordPair::Unpaired), Vec::new())
        };
        return (0..num_pipes)
            .map(|pid| {
                let cmd = if paired {
                    template.replace("{R1}", &a[pid]).replace("{R2}", &b[pid])
                } else {
                    template.replace("{}", &a[pid])
                }
                .replace("{n}", &pid.to_string());
                sh(&cmd)
            })
            .collect();
    }
    if template.contains("{n}") {
        warn!(
            "{{n}} is not expanded by --exec-batch; did you mean -x/--exec (one command per pipe)?"
        );
    }
    let cmd = if paired {
        let (r1s, r2s) = (paths(RecordPair::R1), paths(RecordPair::R2));
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
        template.replace("{}", &paths(RecordPair::Unpaired).join(" "))
    };
    Ok(vec![sh(&cmd)?])
}
