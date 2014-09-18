from subprocess import Popen, PIPE

queueopts = {
        "short6" : "cput=1:00:00,walltime=1:00:00",
        "medium6" : "cput=6:00:00,walltime=6:00:00",
        "long6" : "cput=24:00:00,walltime=24:00:00",
        "vlong6" : "cput=120:00:00,walltime=120:00:00"
        }

def pbssubmit(name, cmd, queue="short6", outfile=None):
    print "submitting job \"%s\" to the batch system:\n%s" % (name, cmd)
    qsub = Popen("qsub", stdin=PIPE, stdout=PIPE, stderr=PIPE)

    lines = ["#PBS -S /bin/bash",   # shell to use
            "#PBS -j oe",           # merge stderr with stdout
            "#PBS -N %s" % name,    # job name
            "#PBS -q %s" % queue,   # queue name
            "#PBS -l %s" % queueopts[queue]]    # cpu and wall time options

    if outfile:
        lines.append("#PBS -o %s" % outfile)

    lines.append("cd ${PBS_O_WORKDIR}")
    lines.append(cmd)

    pbsstring = '\n'.join(lines)

    qsub.stdin.write(pbsstring)
    qsub.stdin.close()

    return qsub
