# Imports
import os
import argparse
import glob
import subprocess

# Main class
class CopyGenotyping(object):
    def __init__(self):
        self.destination = '/srv/gs1/projects/scg/Bina/VA/IL-calls/ILarray'
        self.server = 'scg3'

    # Primary function to be executed
    def copy_files(self, path, server=None, destination=None):
        if not self.check_path(path):
            print "Drive not found. %s" % path
            exit(1)
        if server:
            self.server = server
        if destination:
            self.destination = destination
        genome_dirs = self.list_contents(path)
        self.copy_genotyping(genome_dirs)

    # Show all contents of drive matching LP*
    # Genomes from Illumina follow the naming convention like this: LP6005692-DNA_B12
    def list_contents(self, path):
        genome_glob = "%s/LP*" % path
        genome_dirs = glob.glob(genome_glob)
        return genome_dirs

    # Generate the path to the genotyping file
    def get_genotyping_path(self, sample_barcode):
        genotyping_file = "FinalReport_HumanOmni25M-8v1-1_%s.txt" % sample_barcode
        path = os.path.join('Genotyping', genotyping_file)
        return path

    # Get the path to the genotyping file and run the copy command
    def copy_genotyping(self, genome_list):
        for genome_path in genome_list:
            genome = os.path.basename(genome_path)
            genotyping_path = genome_path + "/" + self.get_genotyping_path(genome)
            if not self.check_path(genotyping_path):
                print "Genotyping file not present! %s" % genotyping_path
                return None
            self.run_rsync(genotyping_path)

    # Copy command
    def run_rsync(self, genotyping_path):
        destination = self.server + ":" + self.destination
        print "Copying %s to %s" % (genotyping_path, destination)
        cmd = ["rsync", "-avLP", genotyping_path, destination]
        return_code = subprocess.call(cmd)#, shell=True)
        if not self.check_return_code(return_code):
            print "Copy of %s failed!" % genotyping_path

    # Check exit code from a job
    def check_return_code(self, return_code):
        if return_code != 0:
            return False
        return True

    # Check file path
    def check_path(self, path):
        if not os.path.exists(path):
            return False
        return True


def parse_command_line():
    parser = argparse.ArgumentParser(
        description = 'This script copies genotyping files from Illumina hard drives to scg3.')

    parser.add_argument("--drive_path", type=str, default=None,
                                help="Path to drive containing files to be copied.")
    parser.add_argument("--server", type=str, default=None,
                                help="Set destination server.  SCG3 default.")
    parser.add_argument("--destination", type=str, default=None,
                                help="Set destination path on the server.  Default mvp genotyping location.")

    options = parser.parse_args()
    if options.drive_path is None:
        print "Drive path is required to proceed.  Use --drive_path to specify a hard drive."
        exit(0)
    return options

if __name__ == "__main__":
    options = parse_command_line()
    copy = CopyGenotyping()
    copy.copy_files(options.drive_path)