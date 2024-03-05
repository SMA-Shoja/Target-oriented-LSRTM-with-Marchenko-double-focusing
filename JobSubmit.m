% Script for submitting Target-oriented LSRTM code to a cluster.
% Adjust it to your cluser specifications.

ip  = java.net.InetAddress.getLocalHost.getHostAddress().string;
pctconfig('hostname',ip);
totalNumberOfWorkers = 60;
myCluster = parcluster('DelftBlue');
myCluster.AdditionalProperties.AdditionalSubmitArgs= '--job-name=TOLSRTM --time=06:00:00 --partition=memory --account=research-ceg-gse --nodes=3 --ntasks-per-node=20 --cpus-per-task=1 --mem=600GB';
job = batch(myCluster,'TargetOrientedLSRTM',2,'Pool',totalNumberOfWorkers-1,'CurrentFolder','.');