#### AWS ecr login
#### Assume aws credentials has been already configured. Please check your ~/.aws/credentials
aws ecr get-login-password --region us-west-2 | docker login --username AWS --password-stdin 204154409870.dkr.ecr.us-west-2.amazonaws.com

#### build docker images
# docker image build -t 204154409870.dkr.ecr.us-west-2.amazonaws.com/bwa:latest -f docker/bwa.dockerfile .
# docker image build -t 204154409870.dkr.ecr.us-west-2.amazonaws.com/cutadapt:latest -f docker/cutadapt.dockerfile .
# docker image build -t 204154409870.dkr.ecr.us-west-2.amazonaws.com/featurecounts:latest -f docker/featurecounts.dockerfile .
# docker image build -t 204154409870.dkr.ecr.us-west-2.amazonaws.com/macs2:latest -f docker/macs2.dockerfile .
# docker image build -t 204154409870.dkr.ecr.us-west-2.amazonaws.com/rseqc:latest -f docker/rseqc.dockerfile .
# docker image build -t 204154409870.dkr.ecr.us-west-2.amazonaws.com/samplot:latest -f docker/samplot.dockerfile .
# docker image build -t 204154409870.dkr.ecr.us-west-2.amazonaws.com/scanalysis:latest -f docker/scanalysis.dockerfile .
# docker image build -t 204154409870.dkr.ecr.us-west-2.amazonaws.com/seurat:latest -f docker/seurat.dockerfile .
# docker image build -t 204154409870.dkr.ecr.us-west-2.amazonaws.com/skbio:latest -f docker/skbio.dockerfile .
# docker image build -t 204154409870.dkr.ecr.us-west-2.amazonaws.com/star:latest -f docker/star.dockerfile .
# docker image build -t 204154409870.dkr.ecr.us-west-2.amazonaws.com/umi_tools:latest -f docker/umi_tools.dockerfile .
# docker image build -t 204154409870.dkr.ecr.us-west-2.amazonaws.com/bamtofrag:latest -f docker/bamtofrag.dockerfile .

#  For a new image, after building the image, log into AWS console > amazon ECR > create repository with a new image name ex,"samplot".
#  Then push the newly built image to ECR as below

####  push images to AWS ECR
# docker push 204154409870.dkr.ecr.us-west-2.amazonaws.com/bwa:latest
# docker push 204154409870.dkr.ecr.us-west-2.amazonaws.com/cutadapt:latest
# docker push 204154409870.dkr.ecr.us-west-2.amazonaws.com/featurecounts:latest
# docker push 204154409870.dkr.ecr.us-west-2.amazonaws.com/macs2:latest
# docker push 204154409870.dkr.ecr.us-west-2.amazonaws.com/rseqc:latest
# docker push 204154409870.dkr.ecr.us-west-2.amazonaws.com/samplot:latest
# docker push 204154409870.dkr.ecr.us-west-2.amazonaws.com/scanalysis:latest
# docker push 204154409870.dkr.ecr.us-west-2.amazonaws.com/seurat:latest
# docker push 204154409870.dkr.ecr.us-west-2.amazonaws.com/skbio:latest
# docker push 204154409870.dkr.ecr.us-west-2.amazonaws.com/star:latest
# docker push 204154409870.dkr.ecr.us-west-2.amazonaws.com/umi_tools:latest
# docker push 204154409870.dkr.ecr.us-west-2.amazonaws.com/bamtofrag:latest


####  pull images from AWS ECR
# docker pull 204154409870.dkr.ecr.us-west-2.amazonaws.com/bwa:latest
# docker pull 204154409870.dkr.ecr.us-west-2.amazonaws.com/cutadapt:latest
# docker pull 204154409870.dkr.ecr.us-west-2.amazonaws.com/featurecounts:latest
# docker pull 204154409870.dkr.ecr.us-west-2.amazonaws.com/macs2:latest
# docker pull 204154409870.dkr.ecr.us-west-2.amazonaws.com/rseqc:latest
# docker pull 204154409870.dkr.ecr.us-west-2.amazonaws.com/samplot:latest
# docker pull 204154409870.dkr.ecr.us-west-2.amazonaws.com/scanalysis:latest
# docker pull 204154409870.dkr.ecr.us-west-2.amazonaws.com/seurat:latest
# docker pull 204154409870.dkr.ecr.us-west-2.amazonaws.com/skbio:latest
# docker pull 204154409870.dkr.ecr.us-west-2.amazonaws.com/star:latest
# docker pull 204154409870.dkr.ecr.us-west-2.amazonaws.com/umi_tools:latest
# docker pull 204154409870.dkr.ecr.us-west-2.amazonaws.com/bamtofrag:latest

