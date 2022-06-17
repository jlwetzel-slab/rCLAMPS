# Script to run HMMER 3 on protein sequences and parse the outputs

import sys, os, re, string, random, operator, gzip
from copy import deepcopy

if sys.platform == 'darwin':
    HOMEDIR = '/Users/jlwetzel/' #for MAC OS X
elif sys.platform == 'linux2':
    HOMEDIR = '/home/jlwetzel/'  #for kernel 2.x or 3.x

#------------------------------------------------------------------------------#
#--------------------- Run HMMER & Parse Raw Output Files ---------------------#
#------------------------------------------------------------------------------#

def runhmmer232(outfile, hmmfile, protfile, hmmname, descs):
    """Runs HMMER 2.3.2 and parses results, storing them into the outfile."""
    #if not hmmfile.endswith('-2'):
    #  os.system('hmmconvert -a '+hmmfile+' '+hmmfile+'-2')
    #  os.system('chmod 755 '+hmmfile+'-2')
    #  hmmfile = hmmfile+'-2'
    os.system("hmmsearch --cut_ga --informat fasta "+hmmfile+" "+protfile+" > "+outfile)
    # os.system('rm '+hmmfile)
    parsehmmer232(outfile, hmmname, descs)


def runhmmer3(outfile, hmmfile, protfile, hmmname, descs, hmmerDir = None):
    """Runs HMMER 3 and parses results, storing them into the outfile."""
    
    #print hmmerDir
    if hmmerDir is None:
      syscall = "hmmsearch --notextw --cut_ga "+hmmfile+" "+protfile+" > "+outfile
    elif hmmerDir is not None:
      syscall = hmmerDir+"./hmmsearch " + \
        "--notextw --cut_ga "+hmmfile+" "+protfile+" > "+outfile
    """
    elif sys.platform == 'darwin':
      syscall = "/Users/jlwetzel/src/hmmer-3.1b2/bin/./hmmsearch " + \
        "--notextw --cut_ga "+hmmfile+" "+protfile+" > "+outfile
    elif sys.platform == 'linux2':
      syscall = HOMEDIR+"/src/hmmer-3.3.1/src/./hmmsearch " + \
        "--notextw --cut_ga "+hmmfile+" "+protfile+" > "+outfile
      print(syscall)
    """
    #syscall = "/home/jlwetzel/src/hmmer-3.1b2/bin/./hmmsearch " + \
    #  "-Z 10 --notextw --domT 0 "+hmmfile+" "+protfile+" > "+outfile
    #syscall = "/usr/local/bin//hmmsearch " + \
    #  "-Z 10 --notextw --domT 0 "+hmmfile+" "+protfile+" > "+outfile

    # [Kaiqian] lrwxr-xr-x  1 kaiqianz  wheel  33 Mar 17 09:39 /usr/local/bin/hmmsearch -> ../Cellar/hmmer/3.3/bin/hmmsearch

    #print syscall
    os.system(syscall)
    #os.system('cp %s %s' %(outfile,outfile+'.orig'))
    parsehmmer3(outfile, hmmname, descs)


def parsehmmer232(filename, hmmname, descs):
  """Parse output from HMMER 2.3.2 into a readable, space-delineated
  table with all important information."""

  tmpfile = '/tmp/' + idgen() + '.txt'
  OUTFILE = open(tmpfile, 'w')

  # Some variables to be used while looping through file
  reach    = False # Have we reached a match yet?
  multihmm = False # Does the HMM match span multiple lines?

  cid = '' # Current ID that we're looking at (for matching purposes)
  for l in open(filename):
    i = l.strip().split()
    if len(i) < 1: continue # No pertinent information here

    # If i[0] is an ID that we're interested in:
    if len(i) > 2 and i[0][:-1] in descs and i[0][-1]==':' and i[1]=='domain':
      assert not reach, 'Could not correctly parse ' + filename
      cid = i[0][:10]

      reach = True # We have reached a match

      # Parse all possible info from this line:
      protID   = i[0][:-1]
      bitscore = i[10][:-1]
      evalue   = i[13]
      alistart = i[6]
      aliend   = i[8][:-1]
      alimatch = ''

    elif reach:
      if i[0][:3] == '*->': # Start of HMM match
        if '<' not in i[0]: # Multi-line HMM match
          multihmm = True
          hmmmatch = i[0][3:]
        else:
          hmmmatch = i[0][3:i[0].find('<')]

      elif i[0] == cid and len(i) > 3: # Part of an alignment match
        if len(alimatch) < 1: alimatch = i[2]
        else:                 alimatch += i[2]

        if aliend == i[3]: # Reached the end
          OUTFILE.write('\t'.join([protID, hmmname, evalue, bitscore, alistart,
                                  aliend, hmmmatch, alimatch, descs[protID] \
                                  if protID in descs else ''])+'\n')
          reach = False
          multihmm = False

      elif multihmm and (len(i)==1 or '+' not in l):
        if '<' in i[0]:
          multihmm = False
          hmmmatch += i[0][:i[0].find('<')]
        elif i[0] != 'CS': # Not actually sure what I needed this for...
          hmmmatch += i[0]
  OUTFILE.close()

  # Overwrite the original HMMER output:
  NOUT = open(filename,'w')
  NOUT.write('\t'.join(['#Target','E-value','BitScore','TargetStart','TargetEnd','HMM_Seq','Ali_Seq','Description'])+'\n')
  with open(tmpfile) as x: map(NOUT.write, sorted(x.readlines()))
  NOUT.close()
  os.system('rm '+tmpfile)


def parsehmmer3(filename, hmmname, descs):
  """Parses output from a HMMER 3 output file into a readable table format
  with all necessary information."""

  tmpfile = '/tmp/' + idgen() + '.txt'
  OUTFILE = open(tmpfile, 'w')

  reach = False # Have we reached a match yet?
  hmmstart,protID = '',''

  #print descs
  for l in open(filename):
    i = l.strip().split()
    if len(i) < 1: continue # No pertinent information here

    if i[0] == '==':
      # We've reached one domain already, so print it out.
      if reach:
        OUTFILE.write('\t'.join([protID,hmmname,evalue,bitscore,hmmstart,hmmend,alistart,aliend,hmmmatch,alimatch,descs[protID] if protID in descs else ''])+'\n')
        hmmstart,protID = '',''

      reach = True
      bitscore = i[4]
      evalue   = i[8]

    elif reach: # Otherwise, parse out information to augment domain
      if i[0].startswith('>>') or i[0].startswith('Internal'):
        OUTFILE.write('\t'.join([protID,hmmname,evalue,bitscore,hmmstart,hmmend,alistart,aliend,hmmmatch,alimatch,descs[protID] if protID in descs else ''])+'\n')
        reach = False
        hmmstart,protID = '',''

      elif i[0].startswith(hmmname):
        if hmmstart == '':
          hmmstart  = i[1]
          hmmmatch  = i[2]
        else:
          hmmmatch += i[2]
        hmmend = i[3]
      elif i[0] in descs:
        if protID == '':
          protID    = i[0] # if '|' not in i[0] else i[0].split('|')[1]
          alistart  = i[1]
          alimatch  = i[2]
        else:
          alimatch += i[2]
        aliend = i[3]
  OUTFILE.close()

  #print tmpfile
  # Overwrite the original HMMER output:
  NOUT = open(filename,'w')
  NOUT.write('\t'.join(['#Target','E-value','BitScore','HMM_Start','HMM_End','TargetStart','TargetEnd','HMM_Seq','Ali_Seq','Description'])+'\n')
  with open(tmpfile) as x: map(NOUT.write, sorted(x.readlines()))
  NOUT.close()
  os.system('rm '+tmpfile)

#------------------------------------------------------------------------------#
#------------------------------ Helper Functions ------------------------------#
#------------------------------------------------------------------------------#

def idgen(size=10, chars=string.ascii_letters + string.digits):
    """Generates a random string (for use as a temporary file) of size n and
    from uppercase, lowercase ASCII letters and numbers."""
    return ''.join(random.choice(chars) for x in range(size))


def getdescs(filename):
  """Given a protein filename, return a dictionary of FASTA descriptions
  (with protein ID as the key)."""

  if filename.endswith('gz'):
    lines = [a for a in gzip.open(filename) if a.startswith('>')]
  else:
    lines = [a for a in open(filename) if a.startswith('>')]

  return {a.strip().split()[0][1:]:' '.join(a.strip().split()[1:]) for a in lines}


def sorttable(table, cols):
    for col in reversed(cols):
        table = sorted(table, key=operator.itemgetter(col))
    return table
