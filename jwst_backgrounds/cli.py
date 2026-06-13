"""
Usage:
    jwst_backgrounds <ra> <dec> <wavelength> [options]
    jwst_backgrounds --field=<name> <wavelength> [options]
    jwst_backgrounds --list-fields

Options:
    -h --help                     show this help
    --field=<name>                use a named JWST field instead of <ra> <dec>
    --list-fields                 list the known field names and exit
    --cache=<ver>                 STScI cache version: 2.0 or 3.0 [default: 2.0]
    --thresh=<float>              threshold factor relative to the minimum background [default: 1.1]
    --day=<integer>               which day in the year for which to extract the background
    --showsubbkgs                 show background components in the bathtub plot
    --background_file=<string>    output file name for the background [default: background.txt]
    --bathtub_file=<string>       output file name for the bathtub curve [default: background_versus_day.txt]

Help:
    For help using this tool, please contact the jwst help desk at jwsthelp.stsci.edu

"""

from jwst_backgrounds import jbt
from jwst_backgrounds.docopt import docopt


def main():
    """Main CLI entrypoint."""
    opt = docopt(__doc__)

    if opt['--list-fields']:
        from jwst_backgrounds import fields
        print("Known JWST fields (name, RA, DEC, cache radius deg):")
        for name in fields.list_fields():
            ra, dec, rad = fields.DEEP_FIELDS[name]
            print("  {0:12s} {1:11.5f} {2:11.5f}  {3:.2f}".format(name, ra, dec, rad))
        return

    if opt['--field']:
        from jwst_backgrounds import fields
        _, ra, dec, _ = fields.resolve(opt['--field'])
    else:
        ra = float(opt['<ra>'])
        dec = float(opt['<dec>'])
    wavelength = float(opt['<wavelength>'])

    thisday = None if opt['--day'] is None else int(opt['--day'])
    cache_url = jbt.CACHE_URLS.get(opt['--cache'], jbt.DEFAULT_CACHE_URL)

    jbt.get_background(ra, dec, wavelength,
                       thresh=float(opt['--thresh']), thisday=thisday, showsubbkgs=opt['--showsubbkgs'],
                       background_file=opt['--background_file'], bathtub_file=opt['--bathtub_file'],
                       cache_url=cache_url)
