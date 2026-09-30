use std::path::Path;

use anyhow::Result;
use binseq::BinseqRecord;
use serde_json::Value;

use crate::commands::qc::{
    base_content::PerBaseSequenceContent, base_quality::PerBaseSequenceQuality,
    dup_levels::SequenceDuplicationLevels, gc_content::PerSequenceGcContent,
    seq_length::SequenceLengthDistribution, seq_quality::PerSequenceQuality,
};

#[derive(Clone)]
pub enum QcModuleType {
    BaseQuality(PerBaseSequenceQuality),
    SeqQuality(PerSequenceQuality),
    BaseContent(PerBaseSequenceContent),
    GcContent(PerSequenceGcContent),
    SeqLength(SequenceLengthDistribution),
    Duplication(SequenceDuplicationLevels),
}
impl QcModuleType {
    pub fn push<R: BinseqRecord>(&mut self, record: &R) {
        match self {
            Self::BaseQuality(x) => x.push(record),
            Self::SeqQuality(x) => x.push(record),
            Self::BaseContent(x) => x.push(record),
            Self::GcContent(x) => x.push(record),
            Self::SeqLength(x) => x.push(record),
            Self::Duplication(x) => x.push(record),
        }
    }
    /// Called once, when a thread has finished all of its batches: merges
    /// thread-local state into shared state.
    pub fn sync_final(&mut self) {
        match self {
            Self::BaseQuality(x) => x.sync_final(),
            Self::SeqQuality(x) => x.sync_final(),
            Self::BaseContent(x) => x.sync_final(),
            Self::GcContent(x) => x.sync_final(),
            Self::SeqLength(x) => x.sync_final(),
            Self::Duplication(x) => x.sync_final(),
        }
    }
    pub fn finish(&mut self, outdir: &Path) -> Result<()> {
        match self {
            Self::BaseQuality(x) => x.finish(outdir),
            Self::SeqQuality(x) => x.finish(outdir),
            Self::BaseContent(x) => x.finish(outdir),
            Self::GcContent(x) => x.finish(outdir),
            Self::SeqLength(x) => x.finish(outdir),
            Self::Duplication(x) => x.finish(outdir),
        }
    }
    /// Renders this module's headline stats as a markdown section.
    pub fn summarize(&self) -> String {
        match self {
            Self::BaseQuality(x) => x.summarize(),
            Self::SeqQuality(x) => x.summarize(),
            Self::BaseContent(x) => x.summarize(),
            Self::GcContent(x) => x.summarize(),
            Self::SeqLength(x) => x.summarize(),
            Self::Duplication(x) => x.summarize(),
        }
    }
    /// This module's JSON key and structured stats (`Null` if it has nothing to report).
    pub fn json(&self) -> (&'static str, Value) {
        match self {
            Self::BaseQuality(x) => ("per_base_quality", x.json()),
            Self::SeqQuality(x) => ("per_sequence_quality", x.json()),
            Self::BaseContent(x) => ("per_base_content", x.json()),
            Self::GcContent(x) => ("per_sequence_gc", x.json()),
            Self::SeqLength(x) => ("sequence_length", x.json()),
            Self::Duplication(x) => ("duplication", x.json()),
        }
    }
}
