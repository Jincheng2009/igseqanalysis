#==============================================================================
#     Copyright (C) 2017  MedImmune, LLC
#     
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#     
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#     
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
#==============================================================================

#igblastn -germline_db_V /rbxpkg/projects/immune_replica/imgt/tcr/human_TRV.fa \
#                           -germline_db_D /rbxpkg/projects/immune_replica/imgt/tcr/human_TRBD.fa \
#               -germline_db_J /rbxpkg/projects/immune_replica/imgt/tcr/human_TRJ.fa \
#               -auxiliary_data /rbxpkg/projects/immune_replica/imgt/optional_file/human_gl.aux \
#               -ig_seqtype TCR \
#               -organism human \
#               -num_alignments_V 1 \
#               -num_alignments_J 1 \
#               -num_alignments_D 1 \
#               -show_translation

import os
import subprocess
import sys

def main():
  entry_path = os.path.dirname(os.path.abspath(__file__))

  vpath = os.path.join(entry_path, "imgt/tcr/human_TRV.fa")
  dpath = os.path.join(entry_path, "imgt/tcr/human_TRBD.fa")
  jpath = os.path.join(entry_path, "imgt/tcr/human_TRJ.fa")
  aux_path = os.path.join(entry_path, "imgt/optional_file/human_gl.aux")

  print(entry_path)

  cmd = ['igblastn', '-germline_db_V', vpath,
                     '-germline_db_D', dpath,
                     '-germline_db_J', jpath,
                     '-ig_seqtype', 'TCR',
                     '-organism', 'human',
                     '-auxiliary_data', aux_path,
                     '-show_translation',
                     '-num_alignments_V', '1',
                     '-num_alignments_D', '1',
                     '-num_alignments_J', '1',]
                     
  ps = subprocess.Popen(cmd, stdout=subprocess.PIPE)

  # Grab stdout line by line as it becomes available.  This will loop until 
  # p terminates.
  while ps.poll() is None:
      l = ps.stdout.readline() # This blocks until it receives a newline.
      sys.stdout.write(l)
  # When the subprocess terminates there might be unconsumed output 
  # that still needs to be processed.
  sys.stdout.write(ps.stdout.read())
  ps.stdout.close()

if __name__ == "__main__":
    main()
