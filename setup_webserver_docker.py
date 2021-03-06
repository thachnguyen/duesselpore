#from netifaces import interfaces, ifaddresses, AF_INET
import os, sys


def ip4_addresses():
    ip_list = []
    for interface in interfaces():
        try:
            for link in ifaddresses(interface)[AF_INET]:
                ip_list.append(link['addr'])
        except Exception:
            pass
    return ip_list

if __name__=="__main__":
	#print('setup IP address')
	# iplist = ip4_addresses()
	# f = open('/home/ag-rossi/projects/duesselpore/NGS_webserver/settings.py', 'r').readlines()
	# f[27] = f[27][:-2]+','+ str(iplist)[1:] + '\n'
	# f1 = open('/home/ag-rossi/projects/duesselpore/NGS_webserver/settings.py', 'w')
	# f1.writelines(f)
	# f1.close()
	print('Updating duesselpore')
	os.system('git -C /home/ag-rossi/projects/duesselpore pull')

	if len(sys.argv) >1:
		if sys.argv[1]=='light':
			print('Downloading human reference genome')
			os.system('wget ftp://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -P /home/ag-rossi/ReferenceData/')
			os.system('wget http://ftp.ensembl.org/pub/release-102/gtf/homo_sapiens/Homo_sapiens.GRCh38.102.gtf.gz -P /home/ag-rossi/ReferenceData/')
			os.system('wget ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -P /home/ag-rossi/ReferenceData/')
		print('creating human reference genome indexes, please wait')
		os.system('minimap2 -t 4 -k14 -w10 -d /home/ag-rossi/ReferenceData/reference_human.mmi /home/ag-rossi/ReferenceData/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz')
		os.unlink('/home/ag-rossi/ReferenceData/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz')
		
		if sys.argv[1]=='full':
			print('Downloading human reference genome')
			os.system('wget ftp://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -P /home/ag-rossi/ReferenceData/')
			os.system('wget http://ftp.ensembl.org/pub/release-102/gtf/homo_sapiens/Homo_sapiens.GRCh38.102.gtf.gz -P /home/ag-rossi/ReferenceData/')
			os.system('wget ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -P /home/ag-rossi/ReferenceData/')
			print('creating human reference genome indexes, please wait')
			os.system('minimap2 -t 4 -k15 -w10 -d /home/ag-rossi/ReferenceData/reference_human.mmi /home/ag-rossi/ReferenceData/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz')

			os.unlink('/home/ag-rossi/ReferenceData/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz')


			print('Downloading mouse reference genome')
			os.system('wget ftp://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz -P /home/ag-rossi/ReferenceData/')
			os.system('wget ftp://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz -P /home/ag-rossi/ReferenceData/')
			os.system('wget http://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz -P /home/ag-rossi/ReferenceData/')
			print('creating mouse reference genome indexes, please wait')
			os.system('minimap2 -t 4 -k15 -w10 -d /home/ag-rossi/ReferenceData/Mus_musculus.mmi /home/ag-rossi/ReferenceData/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz')

			os.unlink('/home/ag-rossi/ReferenceData/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz')
			

			print('Downloading rat reference genome')
			os.system('wget ftp://ftp.ensembl.org/pub/release-102/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz -P /home/ag-rossi/ReferenceData/')
			os.system('wget ftp://ftp.ensembl.org/pub/release-102/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.102.gtf.gz -P /home/ag-rossi/ReferenceData/')
			os.system('wget http://ftp.ensembl.org/pub/release-102/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.Rnor_6.0.cdna.all.fa.gz -P /home/ag-rossi/ReferenceData/')
			print('creating rat reference genome indexes, please wait')
			os.system('minimap2 -t 4 -k15 -w10 -d /home/ag-rossi/ReferenceData/Rattus_norvegicus.mmi /home/ag-rossi/ReferenceData/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz')

			os.unlink('/home/ag-rossi/ReferenceData/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz')

			print('Downloading zebrafish reference genome')
			os.system('wget ftp://ftp.ensembl.org/pub/release-102/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz -P /home/ag-rossi/ReferenceData/')
			os.system('wget ftp://ftp.ensembl.org/pub/release-102/gtf/danio_rerio/Danio_rerio.GRCz11.102.gtf.gz -P /home/ag-rossi/ReferenceData/')
			os.system('wget http://ftp.ensembl.org/pub/release-102/fasta/danio_rerio/cdna/Danio_rerio.GRCz11.cdna.all.fa.gz -P /home/ag-rossi/ReferenceData/')
			print('creating zebrafish reference genome indexes, please wait')
			os.system('minimap2 -t 4 -k15 -w10 -d /home/ag-rossi/ReferenceData/Danio_rerio.mmi /home/ag-rossi/ReferenceData/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz')

			os.unlink('/home/ag-rossi/ReferenceData/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz')

			print('Downloading C elegans reference genome')
			os.system('wget ftp://ftp.ensembl.org/pub/release-102/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz -P /home/ag-rossi/ReferenceData/')
			os.system('wget ftp://ftp.ensembl.org/pub/release-102/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.102.gtf.gz -P /home/ag-rossi/ReferenceData/') 
			os.system('wget http://ftp.ensembl.org/pub/release-104/fasta/caenorhabditis_elegans/cdna/Caenorhabditis_elegans.WBcel235.cdna.all.fa.gz -P /home/ag-rossi/ReferenceData/')
			print('creating Celegans reference genome indexes, please wait')
			os.system('minimap2 -t 4 -k15 -w10 -d /home/ag-rossi/ReferenceData/Caenorhabditis_elegans.mmi /home/ag-rossi/ReferenceData/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz')
			
			print('Downloading Covid19 reference genome')
			os.system('wget http://ftp.ebi.ac.uk/ensemblgenomes/pub/viruses/fasta/sars_cov_2/dna/Sars_cov_2.ASM985889v3.dna.toplevel.fa.gz -P /home/ag-rossi/ReferenceData/')
			os.system('wget http://ftp.ebi.ac.uk/ensemblgenomes/pub/viruses/gtf/sars_cov_2/Sars_cov_2.ASM985889v3.101.gtf.gz -P /home/ag-rossi/ReferenceData/') 
			os.system('wget http://ftp.ebi.ac.uk/ensemblgenomes/pub/viruses/fasta/sars_cov_2/cdna/Sars_cov_2.ASM985889v3.cdna.all.fa.gz -P /home/ag-rossi/ReferenceData/')
			print('creating Covid19 reference genome indexes, please wait')
			os.system('minimap2 -t 4 -k15 -w10 -d /home/ag-rossi/ReferenceData/Covid19.mmi /home/ag-rossi/ReferenceData/Sars_cov_2.ASM985889v3.dna.toplevel.fa.gz')

			os.unlink('/home/ag-rossi/ReferenceData/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz')
			os.unlink('/home/ag-rossi/ReferenceData/Sars_cov_2.ASM985889v3.dna.toplevel.fa.gz')
