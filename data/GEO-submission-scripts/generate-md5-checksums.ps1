# Get-FileHash .\E19\E19_S1_L001_I1_001.fastq.gz -Algorithm MD5
# Get-FileHash .\E19\E19_S1_L001_R1_001.fastq.gz -Algorithm MD5
# Get-FileHash .\E19\E19_S1_L001_R2_001.fastq.gz -Algorithm MD5

# cd C:\Users\JSC_PC\Box\fastq_files_correct\HCKW5DRXX_sample-fastq

Get-ChildItem -Path .\ -Filter *.fastq.gz -Recurse | Get-FileHash -Algorithm MD5