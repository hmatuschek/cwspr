SET(libcwspr_SOURCES
	gen4.f90 jt4.f90 encode4.f90 chkmsg.f90 packjt.f90 getmet4.f90 interleave4.f90
	encode232.f90 entail.f90 hash.f90 hashing.f90
	wqencode.f90 wqdecode.f90 encode232.f90
	deg2grid.f90 grid2deg.f90 azdist.f90 geodist.f90 inter_wspr.f90
	fmtmsg.f90 fix_contest_msg.f90 to_contest_msg.f90
	nhash.c wspr.cc jt4.cc)
SET(libcwspr_HEADERS cwspr.hh
	wspr.hh jt4.hh)

add_library(libcwspr STATIC ${libcwspr_SOURCES})
set_target_properties(libcwspr PROPERTIES OUTPUT_NAME cwspr)
