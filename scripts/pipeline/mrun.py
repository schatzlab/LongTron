#!/usr/bin/env python
#modified from mrun.py to use with long reads (and no assembly)

from __future__ import print_function
import os
import sys
import shutil
import glob
from collections import defaultdict
from operator import itemgetter

SRA_DOWNLOAD_CMD='fastq-dump --split-3 --gzip -O'
#CURL_CMD = 'curl -O -J -L'
CURL_DOWNLOAD_CMD = 'curl -L'
PACBIO_MMAP_SETTINGS = '-ax splice -uf -C5'
TEMPDIR='/home-1/cwilks3@jhu.edu/work3/cwilks/pacbio/all_sra/tmp'

def check_index(prefix):
    for fn in ['%s.%d.ht2' % (prefix, i+1) for i in range(8)]:
        if not os.path.exists(fn):
            raise RuntimeError('Could not find index file "%s"' % fn)


def mkdir_quiet(dr):
    """ Create directories needed to ensure 'dr' exists; no complaining """
    import errno
    if not os.path.isdir(dr):
        try:
            os.makedirs(dr)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


def slurm_header(name, partition, gb, nhours, ntasks_per_node):
    myname = 'mrun' + name
    return """
#!/bin/bash -l
#SBATCH
#SBATCH --job-name={}
#SBATCH --output={}
#SBATCH --error={}
#SBATCH --nodes=1
#SBATCH --partition={}
#SBATCH --mem={}G
#SBATCH --time={}:00:00
#SBATCH --ntasks-per-node={}
""".format(myname, '.' + myname + '.out', '.' + myname + '.err', partition, gb, nhours, ntasks_per_node).strip()

def make_scripts(args, partition='parallel', gb=4, nhours=24, ntasks_per_node=6):
    #make general subdirectory structure for the whole job
    #further subdirectories will be created by individiual tasks under these
    mkdir_quiet(os.path.join(args.output_dir, 'downloads'))
    if 'INDEXES' not in os.environ:
        raise RuntimeError('Set environment variable INDEXES to a directory that has desired indices in hisat2 subdir')
    idx = os.path.join(os.getenv('INDEXES'), 'hisat2', args.idx)
    check_index(idx)

    #custom SRA/ENA processing with manifest from Rail output
    print('#Creating tasks', file=sys.stderr)
    with open(args.accession_file, 'rb') as fh:
        for ln in fh:
            ln = ln.rstrip()
            fields = ln.split('\t')
            #if more than one field, assume a Rail manifest file
            #with name/id in the 5th column
            name = ln
            urls = None
            if len(fields) > 1:
                name = fields[-1]
                urls = fields[0]
                if len(fields) > 3:
                    urls += "," + fields[2]
            if ln == 'accession' or ln == 'uuid':
                continue
            lines = [slurm_header(name, partition, gb, nhours, ntasks_per_node)]
            fn = None
            #TCGA
            if args.tcga_token:
                uuid = name
                lines.append('python mrun.py --name %s --uuid %s --nthreads %d --gb %d --tcga-client \'%s\' --tcga-token %s --output-dir %s --final-dir %s' %
                             (uuid, uuid, ntasks_per_node, gb, args.tcga_client, args.tcga_token, args.output_dir, args.final_dir))
                fn = '%s.sh' % uuid
            #curl (e.g. ENCODE)
            elif len(fields) > 1:
                lines.append('python mrun.py --name %s --nthreads %d --gb %d --output-dir %s --final-dir %s --urls %s' %
                             (name, ntasks_per_node, gb, args.output_dir, args.final_dir, urls))
                fn = '%s.sh' % name
            #SRA
            else:
                lines.append('python mrun.py --name %s --nthreads %d --gb %d --output-dir %s --final-dir %s' %
                             (name, ntasks_per_node, gb, args.output_dir, args.final_dir))
                fn = '%s.sh' % name
            with open(fn, 'wb') as ofh:
                ofh.write('\n'.join(lines) + '\n')
            print('sbatch ' + fn)
                

def run(cmd, desc):
    if not cmd:
        return
    print(cmd, file=sys.stderr)
    ret = os.system(cmd)
    if ret != 0:
        raise RuntimeError(desc + ' failed with exitlevel %d' % ret)

def default_preprocess(args, cwd):
    files1 = sorted(glob.glob(os.path.join(cwd,'*1.fastq*')))
    files2 = sorted(glob.glob(os.path.join(cwd,'*2.fastq*')))
    if len(files1) == 0 and len(files2) == 0:
        files1 = sorted(glob.glob(os.path.join(cwd,'*.fastq*')))
    return([files1, files2, None])

def tcga_preprocess(args, cwd):
    '''Need to decompress/untar downloaded file and determine
        set of paired FASTQ filenames'''
    
    #TCGA has two types (determined by center):
    #1) UNCLCCC type:
    #single non-uuid named .tar.gz file containing 1-2 ?.fastq files
    #where ? is either 1 or 2
    source_file = None
    files = glob.glob(os.path.join(cwd,'*.gz'))
    if len(files) == 1:
        source_file = files[0] 
        run('tar -C %s -zxvf %s' % (cwd, files[0]), 'gunzip_untar')
    #2) BCGSC type:
    #single non-uuid named .tar file containing Nx2 ?.fastq.gz files
    #where ? is either 1 or 2 and N is usually 2
    elif len(files) == 0:
        files = glob.glob(os.path.join(cwd,'*.tar'))
        assert(len(files)==1)
        source_file = files[0] 
        run('tar -C %s -xvf %s' % (cwd, files[0]), 'untar')
   
    (files1, files2, unused) = default_preprocess(args, cwd) 
    return([files1, files2, source_file])
    

def clean_sra_cache(args, sra_cache_dir): 
    sra_temp_fn = os.path.join(sra_cache_dir, args.name, ".sra.cache")
    if os.path.exists(sra_temp_fn):
        os.remove(sra_temp_fn)
    sra_temp_fn = os.path.join(sra_cache_dir, args.name, ".sra.vdbcache.cache")
    if os.path.exists(sra_temp_fn):
        os.remove(sra_temp_fn)


def do_pipeline(args):
    idx = os.path.join('/scratch/groups/blangme2/indexes', 'minimap2', args.idx)
    if os.getenv('INDEXES') is not None:
        idx = os.path.join(os.getenv('INDEXES'), 'minimap2', args.idx)
    download_dir = os.path.join(args.output_dir,'downloads',args.name)
    sam_fn = os.path.join(download_dir,args.name + '.sam')
    bam_fn = os.path.join(download_dir,args.name + '.bam')
    sbam_fn = os.path.join(download_dir,args.name + '.sorted.bam')
    sp_fn = os.path.join(download_dir,args.name + '.splices.tsv.bgz')
    iso_fn = os.path.join(download_dir,args.name + '.raw_isoforms.tsv.gz')
    ver_fn = os.path.join(download_dir,args.name + '.versions')
    temp_dir = os.path.join(download_dir,args.name + '.temp')
    sra_cache_dir = None

    download_cmds = ['%s %s %s' % (SRA_DOWNLOAD_CMD, download_dir, args.name)]
    preprocess = default_preprocess
    #TCGA via GDC
    if args.uuid is not None:
        download_cmds = ['%s download -t %s --log-file %s.dl.log -d %s %s' % \
                    (args.tcga_client, args.tcga_token, \
                    os.path.join(args.output_dir,args.uuid), os.path.join(args.output_dir,'downloads'), args.uuid)]
        preprocess = tcga_preprocess
    #curl (e.g. ENCODE)
    elif args.urls is not None:
        urls = args.urls.split(',')
        #TODO: dont assume these are either compressed or not
        download_cmds = ['%s %s > %s/%s.1.fastq.gz' % (CURL_DOWNLOAD_CMD, urls[0], download_dir, args.name),
                            '%s %s > %s/%s.2.fastq.gz' % (CURL_DOWNLOAD_CMD, urls[1], download_dir, args.name)]
        mkdir_quiet(download_dir)
    #SRA (default)
    else:
        sra_cache_dir = os.path.join(args.output_dir, 'sra')
        mkdir_quiet(download_dir)
    
    #have to do downloading/preprocessing before we can construct
    #later commands due to unknown filenames

    #
    # Download
    #
    try:
        for download_cmd in download_cmds:
            #run(download_cmd, 'download')
            assert os.path.exists(download_dir)
    except Exception as e:
        if not args.keep_intermediates and os.path.exists(download_dir):
            shutil.rmtree(download_dir)
        raise e
    finally:
        if not args.keep_intermediates and sra_cache_dir is not None:
            clean_sra_cache(args, sra_cache_dir)

    #
    # Pre-process
    #
    try:
        fq1s, fq2s, source_file = preprocess(args, download_dir)
    except Exception as e:
        if not args.keep_intermediates and os.path.exists(download_dir):
            shutil.rmtree(download_dir)
        raise e
    if not args.keep_intermediates and source_file is not None:
        os.remove(source_file)


    if len(fq1s) == 0 and args.uuid is None:
        raise RuntimeError('No reads argument')
    if len(fq1s) != len(fq2s) and len(fq1s) > 0 and len(fq2s) > 0:
        raise RuntimeError('Different # of mates in paired FASTQs!')
    for fq in fq1s + fq2s:
        if not os.path.exists(fq):
            raise RuntimeError('No FASTQ file: ' + fq)

    #single mate FASTQ (not paired)
    if len(fq2s) == 0:
        #settings for PacBio
        aligner_cmd = 'minimap2 -t %s %s %s %s \
                    | sambamba view -S -f bam /dev/stdin -o %s' % \
                    (args.nthreads, args.mapper_settings, idx, ','.join(fq1s), bam_fn)
    sb_sort_cmd = 'sambamba sort --tmpdir %s -p -m %dG -t %d -o %s %s' % \
                  (temp_dir, args.gb-1, args.nthreads + args.xthreads, sbam_fn, bam_fn)
    sp_cmd = '/bin/bash -x extract_splices_and_isoforms_from_bam.sh samtools %s %s %s %s %s 2> %s' % \
                     (sbam_fn, sp_fn, iso_fn, download_dir, args.tmp_dir, '%s/extract.err' % download_dir)

    for cmd in [download_cmd, aligner_cmd, sb_sort_cmd, sp_cmd]:
        print(cmd)

    if args.uuid:
        os.system('echo "TCGA download client version" > ' + ver_fn)
        if os.system('%s --version >> %s 2>&1' % (args.tcga_client,ver_fn)) != 0:
            raise RuntimeError('No TCGA download client present')
    os.system('echo "hisat2 version" >> ' + ver_fn)
    if os.system('hisat2 --version >> ' + ver_fn) != 0:
        raise RuntimeError('No hisat2 binary in PATH')
    os.system('echo "sambamba version" >> ' + ver_fn)
    if os.system('sambamba --version 2>> ' + ver_fn) != 0:
        raise RuntimeError('No sambamba binary in PATH')
    os.system('echo "stringtie version" >> ' + ver_fn)
    if os.system('stringtie --version >> ' + ver_fn) != 0:
        raise RuntimeError('No stringtie binary in PATH')


    #
    # ALIGNER + sambamba
    #
    if not args.force and os.path.exists(bam_fn):
        raise RuntimeError('BAM file "%s" exists' % bam_fn)
    try:
        #run(aligner_cmd, 'aligner+sambamba view')
        assert os.path.exists(bam_fn)
    except Exception as e:
        if not args.keep_intermediates and os.path.exists(download_dir):
            shutil.rmtree(download_dir)
        raise e
    if not args.keep_intermediates:
        for fq in fq1s + fq2s:
            os.remove(fq)

    #
    # Sambamba sort
    #
    try:
        mkdir_quiet(temp_dir)
        #run(sb_sort_cmd, 'sambamba sort')
        assert os.path.exists(sbam_fn)
    except Exception as e:
        if not args.keep_intermediates and os.path.exists(download_dir):
            shutil.rmtree(download_dir)
        raise e
    finally:
        if not args.keep_intermediates:
            shutil.rmtree(temp_dir)
    if not args.keep_intermediates:
        os.remove(bam_fn)

    #
    # Extract Splice Jxs
    #
    try:
        run(sp_cmd, 'extract splices + isoforms')
        assert os.path.exists(sp_fn)
        assert os.path.exists(iso_fn)
    except Exception as e:
        if not args.keep_intermediates and os.path.exists(download_dir):
            shutil.rmtree(download_dir)
        raise e
    #if not args.keep_intermediates:
    #    os.remove(sbam_fn)

    run('mv %s %s/' % (sp_fn, args.final_dir), "move_final_file")
    run('mv %s %s/' % (iso_fn, args.final_dir), "move_final_file")
    print('SUCCESS')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Construct pipeline.')
    parser.add_argument('--accession-file', metavar='path', type=str,
			help='file with list of either SRA run accessions (SRR*) or TCGA File UUIDs')
    parser.add_argument('--output-dir', metavar='path', type=str, 
			help='/path/to/root directory where to download and process to; subdirs will be created under this')
    parser.add_argument('--final-dir', metavar='path', type=str, 
			help='/path/to/directory where finished GTF files are moved to')
    parser.add_argument('--tcga-client', metavar='path', type=str, 
            default='/home-1/cwilks3@jhu.edu/gdc-client/python/bin/python /home-1/cwilks3@jhu.edu/gdc-client/bin/gdc-client',
			help='/path/to/Python and GDC client')
    parser.add_argument('--tcga-token', metavar='path', type=str,
			help='/path/to/GDC token file')
    parser.add_argument('--name', metavar='str', type=str, help='Prefix for output files.', default=None)
    parser.add_argument('--uuid', metavar='str', type=str, help='File UUID of TCGA FASTQs', default=None)
    parser.add_argument('--urls', metavar='str', type=str, help='1 or 2 comma separate list of URLs to download from for files other than SRA/TCGA', default=None)
    parser.add_argument('--mapper-settings', metavar='str', type=str, help='mapper/aligner settings', default=PACBIO_MMAP_SETTINGS)
    parser.add_argument('--tmp-dir', metavar='str', type=str, help='temporary directory to use for sorting', default=TEMPDIR)
    parser.add_argument('--gb', metavar='int', type=int, default=100, help='gigabytes of RAM available')
    parser.add_argument('--nthreads', metavar='int', type=int, default=12,
                        help='# threads for hisat')
    parser.add_argument('--xthreads', metavar='int', type=int, default=5,
                        help='additional # threads for sambamba sort/stringtie')
    parser.add_argument('--force', action='store_const', const=True, default=False,
                        help='Overwrite intermediates')
    parser.add_argument('--keep-intermediates', action='store_const', const=True, default=False,
                        help='Don\'t delete intermediate files; if not set, under an exception the whole working directory is removed')
    parser.add_argument('--idx', metavar='str', type=str, help='index prefix', default="hg38")
    args = parser.parse_args()
    if args.accession_file is not None:
        make_scripts(args, gb=args.gb, ntasks_per_node=args.nthreads+args.xthreads)
    elif args.name is not None or args.uuid is not None:
        do_pipeline(args)
    else:
        parser.print_help()
        parser.exit()    
