import sys,os
from sastbx import data_reduction
from sastbx import interface
import iotbx.phil
from iotbx import ranges
from libtbx.utils import Sorry



master_params = iotbx.phil.parse("""\
    tdfXS{

      images {
         background{
           base = None
          .help="Image base name"
          .type=str
          range=None
          .type=str
          .help="Range of images"
          file_name_settings
          .expert_level=1000
          {
            wildcard = "#"
            .help="Serial id wildcard"
            .type=str
            has_leading_zeros=True
            .type=bool
            .help="Does the serial id have leading zeros?"
          }
         }


         dark{
          base = None
          .help="Image base name"
          .type=str
          range=None
          .type=str
          .help="Range of images"
          file_name_settings
          .expert_level=1000
          {
            wildcard = "#"
            .help="Serial id wildcard"
            .type=str
            has_leading_zeros=True
            .type=bool
            .help="Does the serial id have leading zeros?"
          }
         }
         sample{
          base = None
          .help="Image base name"
          .type=str
          range=None
          .type=str
          .help="Range of images"
          file_name_settings
          .expert_level=1000
          {
            wildcard = "#"
            .help="Serial id wildcard"
            .type=str
            has_leading_zeros=True
            .type=bool
            .help="Does the serial id have leading zeros?"
          }
         }
      }


            mask{

        shadow
        .help="Shadow area for offset corrections"
        {
          circle
          .help="A circular shadow area"
          {
            origin = None
            .type= str
            .help="Coordinates of shadow circle"
            radius=None
            .type=float
            .help="Radius of shadow circle"
          }
          polygon
          .multiple=True
          .help="Polygon shadow area"
          {
            corner = None
            .type = str
            .help="Corner of shadow zone. We need at least 3 of these coordinates."
            .multiple = True
          }
        }

        exclude
        .help="Shadow area for offset corrections"
        {
          circle
          .help="A circular shadow area"
          {
            origin = None
            .type= str
            .help="Coordinates of shadow circle"
            radius=None
            .type=float
            .help="Radius of shadow circle"
          }
          polygon
          .multiple=True
          .help="Polygon shadow area"
          {
            corner = None
            .type = str
            .help="Corner of shadow zone. We need at least 3 of these coordinates."
            .multiple = True
          }
        }
      }




     }""")


def get_file_names(base_name, wildcard, has_leading_zeros, user_range_txt):
  file_name_engine = ranges.serial_file_name_handler(base_name, wildcard)
  id_list = ranges.range_to_list( ranges.range_parser(user_range_txt) )
  result = file_name_engine.names_from_range_list( id_list )
  return result




banner = """\n
***************************************************************
*                                                             *
*                                                             *
*                 Processing of 2D fXS data                   *
*                                                             *
*                                                             *
***************************************************************
"""

def helpf(out=None):
  if out is None:
    out = sys.stdout
  print >> out, "NO HELP YET"





def run(args):

  out = sys.stdout

  params = interface.get_input( args, master_params, "tdfXS", banner, helpf)

  bg_list = []# None
  dk_list = []# None
  sm_list = []# None

  if params.tdfXS.images.background.base is not None:
    bg_base  = params.tdfXS.images.background.base
    bg_range = params.tdfXS.images.background.range
    bg_wc    = params.tdfXS.images.background.file_name_settings.wildcard
    bg_hlz   = params.tdfXS.images.background.file_name_settings.has_leading_zeros
    if bg_range is None:
      raise Sorry("Background range cannot be empty")
    bg_list  = get_file_names(bg_base, bg_wc, bg_hlz, bg_range)

  if params.tdfXS.images.dark.base is not None:
    dk_base  = params.tdfXS.images.dark.base
    dk_range = params.tdfXS.images.dark.range
    dk_wc    = params.tdfXS.images.dark.file_name_settings.wildcard
    dk_hlz   = params.tdfXS.images.dark.file_name_settings.has_leading_zeros
    if dk_range is None:
      raise Sorry("Dark range cannot be empty")
    dk_list  = get_file_names(dk_base, dk_wc, dk_hlz, dk_range)

  if params.tdfXS.images.sample.base is not None:
    sm_base  = params.tdfXS.images.sample.base
    sm_range = params.tdfXS.images.sample.range
    sm_wc    = params.tdfXS.images.sample.file_name_settings.wildcard
    sm_hlz   = params.tdfXS.images.sample.file_name_settings.has_leading_zeros
    if sm_range is None:
      raise Sorry("Sample range cannot be empty")

    sm_list  = get_file_names(sm_base, sm_wc, sm_hlz, sm_range)

  print >> out
  print >> out, "----- Processing images -----"
  print >> out
  print >> out, "   Sample Images     : "
  for ii in  sm_list:
    print  >> out, "         ", ii
  print

  print "   Dark Images       : "
  for ii in dk_list:
    print  >> out, "         ",ii
  print

  print "   Background Images : "
  for ii in bg_list:
    print  >> out, "         ",ii
  print





if __name__ == "__main__":
  run(sys.argv[1:])
