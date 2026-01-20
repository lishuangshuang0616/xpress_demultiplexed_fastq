import os
import gzip

os.makedirs("test_input", exist_ok=True)

# Create samplesheet
with open("samplesheet.csv", "w") as f:
    f.write("barcode_id,sequence\n")
    f.write("10,AAAA\n")
    f.write("11,GGGG\n")

# Helper to write fastq
def write_fastq(path, name, seq, qual, compress=False):
    content = f"@{name}\n{seq}\n+\n{qual}\n"
    if compress:
        with gzip.open(path, "wt") as f:
            f.write(content)
    else:
        with open(path, "w") as f:
            f.write(content)

# ID 10 (Compressed)
write_fastq("test_input/test_L01_10_1.fq.gz", "read1_10", "ACGT", "FFFF", compress=True)
write_fastq("test_input/test_L01_10_2.fq.gz", "read1_10", "TGCA", "FFFF", compress=True)

# ID 11 (Uncompressed)
write_fastq("test_input/test_L01_11_1.fq", "read2_11", "TATA", "IIII", compress=False)
write_fastq("test_input/test_L01_11_2.fq", "read2_11", "GCGC", "IIII", compress=False)

print("Test data created.")
