#!/usr/bin/env python
from subprocess import call, PIPE, Popen

def get_tip_rev(text):
    for line in text.splitlines():
        if line.startswith('changeset'):
            rev_num, hashkey = line.split()[1].split(':')
            return rev_num, hashkey

output = Popen(["tests/all_tests.py"], stdout=PIPE, stderr=PIPE, shell=True)

stdout, stderr = output.communicate()
returncode = output.returncode

if returncode != 0:
    print stderr
    result = Popen('hg tip', shell=True, stdout=PIPE)
    stdout, stderr = result.communicate()
    rev_num, hashkey = get_tip_rev(stdout)
    result = Popen('hg strip %s -n' % rev_num, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = result.communicate()
    print 'BAD CHANGESET "%s:%s", FIX BEFORE RESUBMITTING' % (rev_num, hashkey)
    if result.returncode != 0:
        print 'ERROR stripping changset: "%s"\nNOTIFY GAVIN' % stderr
else:
    print stdout

exit(returncode)