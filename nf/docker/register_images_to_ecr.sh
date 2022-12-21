// AWS ecr login
// Assume aws credentials has been already configured. Please check your ~/.aws/credentials
aws ecr get-login-password --region us-west-2 | docker login --username AWS --password-stdin 204154409870.dkr.ecr.us-west-2.amazonaws.com

// build docker images
docker image build -t 204154409870.dkr.ecr.us-west-2.amazonaws.com/bwa:latest - < bwa.dockerfile
docker image build -t 204154409870.dkr.ecr.us-west-2.amazonaws.com/cutadapt:latest - < cutadapt.dockerfile
docker image build -t 204154409870.dkr.ecr.us-west-2.amazonaws.com/featurecounts:latest - < featurecounts.dockerfile
docker image build -t 204154409870.dkr.ecr.us-west-2.amazonaws.com/macs2:latest - < macs2.dockerfile
docker image build -t 204154409870.dkr.ecr.us-west-2.amazonaws.com/rseqc:latest - < rseqc.dockerfile
docker image build -t 204154409870.dkr.ecr.us-west-2.amazonaws.com/scanalysis:latest - < scanalysis.dockerfile
docker image build -t 204154409870.dkr.ecr.us-west-2.amazonaws.com/skbio:latest - < skbio.dockerfile
docker image build -t 204154409870.dkr.ecr.us-west-2.amazonaws.com/star:latest - < star.dockerfile
docker image build -t 204154409870.dkr.ecr.us-west-2.amazonaws.com/umi_tools:latest - < umi_tools.dockerfile
docker image build -t 204154409870.dkr.ecr.us-west-2.amazonaws.com/samplot:latest - < samplot.dockerfile

// For a new image, after building the image, log into AWS console > amazon ECR > create repository with a new image name ex,"samplot".
// Then push the newly built image to ECR as below

// push images to AWS ECR
docker push 204154409870.dkr.ecr.us-west-2.amazonaws.com/bwa:latest
docker push 204154409870.dkr.ecr.us-west-2.amazonaws.com/cutadapt:latest
docker push 204154409870.dkr.ecr.us-west-2.amazonaws.com/featurecounts:latest
docker push 204154409870.dkr.ecr.us-west-2.amazonaws.com/macs2:latest
docker push 204154409870.dkr.ecr.us-west-2.amazonaws.com/rseqc:latest
docker push 204154409870.dkr.ecr.us-west-2.amazonaws.com/scanalysis:latest
docker push 204154409870.dkr.ecr.us-west-2.amazonaws.com/skbio:latest
docker push 204154409870.dkr.ecr.us-west-2.amazonaws.com/star:latest
docker push 204154409870.dkr.ecr.us-west-2.amazonaws.com/umi_tools:latest
docker push 204154409870.dkr.ecr.us-west-2.amazonaws.com/samplot:latest
