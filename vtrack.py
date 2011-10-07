# -*- coding: utf-8 -*-
"""
Provides minimal version tracking header. The fields of the
vheader are:
  call: the parsed call to the script.
  platform: the OS, name of the machine and system version.
  directory: the current directory.
  user: the logged-in user running the script.
  date: date and time of the call.
  'script' SHA1: SHA1 digest of the script as a text file.
  git commit: SHA1 digest of the git commit of the script.

If the script is not properly git-versioned, it is written
in full in the header.

The vheader can be generated by vheader(*sys.argv) anywhere
in the script.

The class vskip is a utility to read-in files with a vheader.
The command 'vskip(open(fname))' will seamlessly strip off
the vheader.

The module exports only two functions, so you can declare
your import as 'from vtrack import vskip, vheader'.
"""

import os
import re
import sys
import datetime
from hashlib import sha1


__author__ = 'Guillaume Filion'
__email__ = 'guillaume.filion@gmail.com'
__version__ = '0.1'



def vheader(script_name, *args, **kwargs):
   """Return a header with version tracking information.
   Pass 'sys.argv' as parameters and set the 'comment_char'
   argument to change the prefix of the header (default '#')."""
   comment_char = kwargs.get('comment_char', '#')
   head_list = []
   # Gather base info.
   for get_base_info in [
         __get_call,
         __get_uname_info,
         __get_directory,
         __get_user_info,
         __get_date,
         __get_python_info
      ]:
      head_list.extend(get_base_info())

   # Gather SHA1 digests.
   for file_name in [script_name] + list(args):
      # Argument files may have a relative or abolute path.
      file_name = file_name if os.path.isabs(file_name) else \
            os.path.join(os.getcwd(), file_name)
      try:
         head_list.extend(__get_sha1(file_name))
      except IOError:
         # Can't open the file (the argument may not be a
         # file name). Just pass.
         pass

   # Gather script git version info.
   try:
      head_list.extend(__get_git_info(script_name))
      append_script = False
   except Exception, e:
      head_list.extend([('git error', str(e))])
      append_script = True

   # Generate the header: concatenate key/value pairs.
   header = ''
   for pair in head_list:
      header += comment_char + ' %s: %s\n' % pair
   if append_script:
      # Append the full script to the header.
      header +=  __get_script(script_name, comment_char) + '\n'

   return header



##  Private functions all return key/value pairs,
##  except __get_script().

def __get_call():
   """Return the system invocation of the script."""
   return [
      ('call', ' '.join(sys.argv))
   ]

def __get_uname_info():
   """Return the name of the machine and the system version."""
   uname = os.uname()
   return [
      ('platform', ' '.join(uname[:3]))
   ]

def __get_directory():
   """Return the path where the script is run."""
   return [
      ('directory', os.getcwd())
   ]

def __get_user_info():
   """Return the user logged in on the terminal."""
   return [
      ('user', os.getlogin())
   ]

def __get_date():
   """Return the date the script is run."""
   return [
      ('date', str(datetime.datetime.today()))
   ]

def __get_python_info():
   """Return the running Python version."""
   return [
      ('python version', '.'.join(str(d) for d in sys.version_info))
   ]

def __get_sha1(file_name):
   """Return the short SHA1 digest of the specified file."""
   sha1_digest = sha1(open(file_name).read()).hexdigest()
   return [
      ('%s SHA1' % file_name, sha1_digest)
   ]

def __get_git_info(script_name):
   """Get git commit. If the file has been changed, or information
   is not available, the whole script is returned."""
   # Will raise ImportError if module git is not installed.
   import git
   # Get repository data. Will raise InvalidGitRepositoryError
   # if the directory is not a repository.
   repo = git.Repo(sys.path[0])
   # Check whether the script is in untracked files.
   if os.path.basename(script_name) in repo.untracked_files:
      raise Exception('%s not tracked' % script_name)
   # Check whether the script is in a committed state.
   changed = [diff.a_blob.name for diff in repo.index.diff(None)]
   if os.path.basename(script_name) in changed:
      raise Exception('%s has been changed since last commit' % script_name)
   # Get the commit's SHA1 digest.
   return [
      ('git commit', repo.head.commit.hexsha)
   ]

def __get_script(script_name, comment_char):
   """Return the full script, prepend the comment character in
   front of every line."""
   comment_char += ' '
   outcommented_script = comment_char + \
      re.sub('\n', '\n' + comment_char, open(script_name).read())
   return outcommented_script



class vskip:
   """File encapsulation class for seamless reading.
   Redefine file iterator to skip lines starting with
   given comment character. Use as 'vskip(open(fname))'."""

   def __init__ (self, finstance, comment_char='#'):
      # Has-a file.
      self.finstance = finstance
      self.comment_char = comment_char

   def __iter__(self):
      """Overwrite the iterator."""
      return self.iterskip()

   def iterskip(self):
      for line in self.finstance:
         if line.startswith(self.comment_char):
            continue
         yield line

   def read(self, *args, **kwargs):
      while True:
         line = self.finstance.readline()
         if not line.startswith(self.comment_char):
            return line + self.finstance.read(*args, **kwargs)

   def readline(self, *args, **kwargs):
      while True:
         line = self.finstance.readline(*args, **kwargs)
         if not line.startswith(self.comment_char):
            return line

   def readlines(self, *args, **kwargs):
      return [
            line \
            for line in self.finstance.readlines(*args, **kwargs) \
            if not line.startswith(self.comment_char)
         ]

   # Encapsulate the file.
   def __getattr__(self, attr):
      return getattr(self.finstance, attr)
