#!/usr/bin/env python

from integrationtest import IntegrationTest, main
import unittest

class TestCreate(IntegrationTest, unittest.TestCase):
    def test_junctions_create(self):
        bam1 = self.inputFiles("test.rnaseq.bam")[0]
        output_file = self.tempFile("create.out")
        print "BAM1 is ", bam1
        for anchor in ["", "30"]:
            expected_file = self.inputFiles("junctions-create/expected-a" +
                                            anchor + ".out")[0]
            if anchor != "":
                anchor = "-a " + anchor
            params = [ "junctions", "create", bam1, anchor, "-o", output_file]
            rv, err = self.execute(params)
            self.assertEqual(rv, 0)
            #self.assertEqual('', err)
            self.assertFilesEqual(expected_file, output_file)

if __name__ == "__main__":
    main()
