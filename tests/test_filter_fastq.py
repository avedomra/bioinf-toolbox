"""
Tests for bioinf_toolbox.filter_fastq and related helpers
"""

import io
import os
import logging
import textwrap
import pytest
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Make sure the project root is importable when running from any directory
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from bioinf_toolbox import filter_fastq, _parse_bounds, _passes_filters, setup_logger


# Helpers

def _make_record(seq: str, qualities: list[int], name: str = "read1") -> SeqRecord:
    """Create a minimal SeqRecord with Phred quality annotations"""
    record = SeqRecord(Seq(seq), id=name, description="")
    record.letter_annotations["phred_quality"] = qualities
    return record


def _write_fastq(path: str, records: list[SeqRecord]) -> None:
    """Write a list of SeqRecords to a FASTQ file"""
    with open(path, "w") as fh:
        SeqIO.write(records, fh, "fastq")


def _read_fastq(path: str) -> list[SeqRecord]:
    """Read all records from a FASTQ file"""
    return list(SeqIO.parse(path, "fastq"))


# Class 1: _parse_bounds

class TestParseBounds:
    def test_single_number_uses_default_lower(self):
        """A single number should set the upper bound; lower stays at default"""
        lo, hi = _parse_bounds(50, default_lower=0.0, default_upper=100.0)
        assert lo == 0.0
        assert hi == 50.0

    def test_tuple_two_values(self):
        """A two-element tuple should be returned as-is"""
        lo, hi = _parse_bounds((20, 80))
        assert lo == 20.0
        assert hi == 80.0


# Class 2: _passes_filters

class TestPassesFilters:
    def test_passes_all_default_filters(self):
        """A typical read should pass with the most permissive defaults"""
        record = _make_record("ATGC", [30, 30, 30, 30])
        assert _passes_filters(record, 0, 100, 0, 2**32, 0) is True

    def test_fails_gc_filter(self):
        """A read with 100 % GC should fail when upper GC bound is 50"""
        record = _make_record("GGGG", [30, 30, 30, 30])
        assert _passes_filters(record, 0, 50, 0, 2**32, 0) is False

    def test_fails_length_filter(self):
        """A read shorter than len_lower should be rejected"""
        record = _make_record("AT", [30, 30])
        assert _passes_filters(record, 0, 100, 10, 2**32, 0) is False

    def test_fails_quality_filter(self):
        """A read with mean quality below threshold should be rejected"""
        record = _make_record("ATGC", [5, 5, 5, 5])
        assert _passes_filters(record, 0, 100, 0, 2**32, 20) is False


# Class 3: filter_fastq – file I/O and integration

class TestFilterFastq:
    def test_output_file_is_created(self, tmp_path):
        """filter_fastq must create the output file inside 'filtered/'"""
        records = [_make_record("ATGCATGC", [30] * 8, "r1")]
        input_file = tmp_path / "input.fastq"
        _write_fastq(str(input_file), records)

        os.chdir(tmp_path)
        filter_fastq(str(input_file), "out.fastq")

        assert (tmp_path / "filtered" / "out.fastq").exists()

    def test_passing_reads_are_written_correctly(self, tmp_path):
        """Reads that pass all filters must appear verbatim in the output"""
        records = [
            _make_record("ATGCATGC", [30] * 8, "pass_read"),
            _make_record("GGGGGGGG", [30] * 8, "fail_gc"),   # 100 % GC, will fail
        ]
        input_file = tmp_path / "input.fastq"
        _write_fastq(str(input_file), records)

        os.chdir(tmp_path)
        filter_fastq(str(input_file), "out.fastq", gc_bounds=(0, 60))

        result = _read_fastq(str(tmp_path / "filtered" / "out.fastq"))
        assert len(result) == 1
        assert result[0].id == "pass_read"

    def test_raises_file_exists_error(self, tmp_path):
        """filter_fastq must raise FileExistsError when output already exists"""
        records = [_make_record("ATGC", [30] * 4)]
        input_file = tmp_path / "input.fastq"
        _write_fastq(str(input_file), records)

        os.chdir(tmp_path)
        filter_fastq(str(input_file), "out.fastq")

        with pytest.raises(FileExistsError):
            filter_fastq(str(input_file), "out.fastq")


# Class 4: logging

class TestLogging:
    def test_log_file_is_created(self, tmp_path):
        """setup_logger must create the log file on the filesystem"""
        log_path = str(tmp_path / "test.log")
        # Remove any cached handlers so a fresh file handler is added
        lg = logging.getLogger("bioinf_toolbox")
        lg.handlers.clear()

        setup_logger(log_path)
        lg.info("test message")

        # Flush handlers
        for h in lg.handlers:
            h.flush()

        assert os.path.exists(log_path)

    def test_error_is_logged_on_existing_output(self, tmp_path, caplog):
        """filter_fastq must log an ERROR when the output file already exists"""
        records = [_make_record("ATGC", [30] * 4)]
        input_file = tmp_path / "input.fastq"
        _write_fastq(str(input_file), records)

        os.chdir(tmp_path)
        filter_fastq(str(input_file), "out.fastq")

        with caplog.at_level(logging.ERROR, logger="bioinf_toolbox"):
            with pytest.raises(FileExistsError):
                filter_fastq(str(input_file), "out.fastq")

        assert any("already exists" in msg for msg in caplog.messages)
