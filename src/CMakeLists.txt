set(cwspr_SOURCES cwspr.cc
	application.cc mainwindow.cc keyer.cc encoder.cc
	waterfall.cc)
set(cwspr_MOC_HEADERS
	application.hh mainwindow.hh keyer.hh encoder.hh
	waterfall.hh)

qt5_wrap_cpp(cwspr_MOC_SOURCES ${cwspr_MOC_HEADERS})
qt5_add_resources(cwsrp_RCC_SOURCES ../shared/resources.qrc)

add_executable(cwspr ${cwspr_SOURCES} ${cwspr_MOC_SOURCES} ${cwsrp_RCC_SOURCES})
target_link_libraries(cwspr libcwspr ${LIBS})
