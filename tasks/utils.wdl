version development

task GetGeneChrPairs {

    input {
        File featureVariantFile
        Array[String]+ columnsToSelect
    }

	command <<<

cat << EOF > script.py
import csv

columns_to_select = '~{sep(",", columnsToSelect)}'.split(",")
with open("~{featureVariantFile}") as f, open("outputPairs.tsv", "w+") as w:
    # probably a bad way to split to get the headers
    headers = f.readline().split(",")
    # this will fail if column isn't in the header
    indices_to_select = [headers.index(h) for h in columns_to_select]

    cut_lines = [[l[i] for i in indices_to_select] for l in csv.reader(f)]
    unique_lines = set('\t'.join(l) for l in cut_lines)
    if unique_lines:
    # this should leave an empty fails
        w.write("\n".join(sorted(unique_lines)))

EOF

python script.py

	>>>

    runtime {
        memory: "2G"
		container: "python:3.10"
    }

    output {
        # [["GeneName1", "chr1"], ["GeneName2","chr2"]]
        Array[Array[String]] output_pairs = read_tsv("outputPairs.tsv")
    }
}
