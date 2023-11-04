"""
	Find closest, previous, next nest dates for a given year, month, day
	Given date could be: year, month, month date or
	year and year day jday= ...
	returned: list for 3 dates: closest, previous, next
"""
import datetime
from datetime import date
from datetime import datetime as dtime
from datetime import timedelta
import mod_time_days
from mod_time_days import get_dtime

def find_nest_date(yr,mo=0,mday=0,dt=4,jday=None):
# HYCOM calendar starts:
	dnst0 = dtime(1900,12,31)
# Define requested date
	if jday is not None:
		dj1 = dtime(yr,1,1)
		dnmb = dj1+timedelta(days=jday-1)
	  yr   = int(dnmb.strftime('%G'))
	  mo   = int(dnmb.strftime('%m'))
	  mday = int(dnmb.strftime('%d'))
	else:
		dnmb = dtime(yr,mo,mday)

	dlt = dnmb-dnst0
	ndays = dlt.days
	d0 = 0.
	while d0 < ndays:
		d0 += dt

	day_nest0 = dnst0+timedelta(days=d0)
	dnest_prev = day_nest0-timedelta(days=dt)
	dnest_next = day_nest0+timedelta(days=dt)

	nest_dates = []
	nest_dates.append([])
	yr_nest   = int(day_nest0.strftime('%G'))
	mo_nest   = int(day_nest0.strftime('%m'))
	md_nest   = int(day_nest0.strftime('%d'))
	jday_nest = int(day_nest0.strftime('%j'))
	dS, dE, dC = get_dtime(yr=yr_nest,mo=mo_nest,mday=md_nest)  # HYCOM date numbers
	nest_dates[0].append((yr_nest,mo_nest,md_nest,jday_nest,dC))

	nest_dates.append([])
	yr_nest   = int(dnest_prev.strftime('%G'))
	mo_nest   = int(dnest_prev.strftime('%m'))
	md_nest   = int(dnest_prev.strftime('%d'))
	jday_nest = int(dnest_prev.strftime('%j'))
	dS, dE, dC = get_dtime(yr=yr_nest,mo=mo_nest,mday=md_nest)
	nest_dates[1].append((yr_nest,mo_nest,md_nest,jday_nest,dC))

	nest_dates.append([])
	yr_nest   = int(dnest_next.strftime('%G'))
	mo_nest   = int(dnest_next.strftime('%m'))
	md_nest   = int(dnest_next.strftime('%d'))
	jday_nest = int(dnest_next.strftime('%j'))
	dS, dE, dC = get_dtime(yr=yr_nest,mo=mo_nest,mday=md_nest)
	nest_dates[2].append((yr_nest,mo_nest,md_nest,jday_nest,dC))


	return 	nest_dates
	




