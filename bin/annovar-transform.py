import argparse
import os

START = 1
END = 2

class AnnotationConversion(object):
    def __init__(self, input, one_to_zero=False):
        # open file
        with open(input) as f:
            for line in f:
                # remove newline
                line = line.rstrip()
                # remove existing commas
                line = self.remove_commas(line)
                # tsv > csv
                line = self.tsv_to_csv(line)
                # Remove null values
                line = self.remove_nulls(line)
                # one to zero conversion
                if one_to_zero is True:
                    line = self.one_to_zero(line)
                # print line
                print line

    def remove_commas(self, line):
        return line.replace(',', ';')

    def tsv_to_csv(self, line):
        return line.replace("\t", ',')

    def remove_nulls(self, line):
        fields = line.split(",")
        for i in range(len(fields)):
            if fields[i] == '.':
                fields[i] = ''
	return ",".join(fields)

    def one_to_zero(self, line):
        fields = line.split(",")
        fields[START] = self.minus_one(fields[START])
        return ",".join(fields)

    def minus_one(self, field):
        try:
            position = str(int(field) - 1)
            return position
        except:
            return field

def parse_command_line():
    parser = argparse.ArgumentParser(
        description = 'This script removes ')

    parser.add_argument("--input", default=None,
                                help="Annovar annotation file")
    parser.add_argument("--one_to_zero", action='store_true', default=False,
                                help="Convert from 1 based to 0 based genomic positions")


    options = parser.parse_args()
    if options.input is None:
        print "Exiting, provide an input file."
        exit(1)

    if not os.path.exists(options.input):
        print "Input file not found. %s" % options.input
        exit(1)


    return options

if __name__ == "__main__":
    options = parse_command_line()
    AnnotationConversion(input=options.input, one_to_zero=options.one_to_zero)
