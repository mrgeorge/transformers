#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import pcorr

pcorr.main(csOption='cen',limOption='bins')
pcorr.main(csOption='sat',limOption='bins')
pcorr.main(csOption='censat',limOption='bins')
pcorr.main(csOption='cen',limOption='complete')
pcorr.main(csOption='sat',limOption='complete')
pcorr.main(csOption='censat',limOption='complete')
#pcorr.main(csOption='cen',limOption='all')
#pcorr.main(csOption='sat',limOption='all')
#pcorr.main(csOption='censat',limOption='all')
#pcorr.main(csOption='cen',limOption='thresh')
#pcorr.main(csOption='sat',limOption='thresh')
#pcorr.main(csOption='censat',limOption='thresh')
