import sys
sys.path.extend(['..', '../src'])

from chippy.prep.mapped_files import mapped_file_handle, internal_instantiation,\
    TEST_DIR_NAME, TEST_NAME, test_instantiation
from cogent.util.unit_test import TestCase, main

class MappedFilesTest(TestCase):
    """we need to test that MappedFiles can be instantiated in another module
    and then work properly here"""
    def test_mapped_files_match_known(self):

        # create from within the Singleton with "inputFile1.fastq"
        mapped_files = test_instantiation()

        # So we know instantiation worked correctly
        test_fastq = '%s.fastq' % TEST_NAME
        self.assertEqual(mapped_files.fastq_fn, test_fastq)

        test_working_dir = '%s-working' % TEST_DIR_NAME
        self.assertEqual(mapped_files.working_dn, test_working_dir)

        full_test_name = '%s/%s' % (test_working_dir, TEST_NAME)
        test_fasta = '%s.fasta' % full_test_name
        self.assertEqual(mapped_files.fasta_fn, test_fasta)

        test_psl = '%s.psl' % full_test_name
        self.assertEqual(mapped_files.psl_fn, test_psl)

        test_contaminated = '%s_contaminated.fastq' % full_test_name
        self.assertEqual(mapped_files.contaminated_fn, test_contaminated)
        test_pristine = '%s_pristine.fastq' % full_test_name
        self.assertEqual(mapped_files.pristine_fn, test_pristine)
        test_combined = '%s_combined.fastq' % full_test_name
        self.assertEqual(mapped_files.combined_fn, test_combined)

        test_combined_sai = '%s_combined.sai' % full_test_name
        self.assertEqual(mapped_files.combined_sai_fn, test_combined_sai)
        test_combined_sam = '%s_combined.sam' % full_test_name
        self.assertEqual(mapped_files.combined_sam_fn, test_combined_sam)

        self.assertEqual(mapped_files.run_record_fn, "%s/run_record.txt" % TEST_DIR_NAME)

        # now try creating a new instance locally with different internals
        mapped_files = mapped_file_handle("inputFile2.fastq", "/blah")

        # Prove we got back the original instance instead
        test_fastq = '%s.fastq' % TEST_NAME
        self.assertEqual(mapped_files.fastq_fn, test_fastq)

        test_working_dir = '%s-working' % TEST_DIR_NAME
        self.assertEqual(mapped_files.working_dn, test_working_dir)

        full_test_name = '%s/%s' % (test_working_dir, TEST_NAME)
        test_fasta = '%s.fasta' % full_test_name
        self.assertEqual(mapped_files.fasta_fn, test_fasta)

        test_psl = '%s.psl' % full_test_name
        self.assertEqual(mapped_files.psl_fn, test_psl)

        test_contaminated = '%s_contaminated.fastq' % full_test_name
        self.assertEqual(mapped_files.contaminated_fn, test_contaminated)
        test_pristine = '%s_pristine.fastq' % full_test_name
        self.assertEqual(mapped_files.pristine_fn, test_pristine)
        test_combined = '%s_combined.fastq' % full_test_name
        self.assertEqual(mapped_files.combined_fn, test_combined)

        test_combined_sai = '%s_combined.sai' % full_test_name
        self.assertEqual(mapped_files.combined_sai_fn, test_combined_sai)
        test_combined_sam = '%s_combined.sam' % full_test_name
        self.assertEqual(mapped_files.combined_sam_fn, test_combined_sam)

        self.assertEqual(mapped_files.run_record_fn, "%s/run_record.txt" % TEST_DIR_NAME)


if __name__ == '__main__':
    main()
    