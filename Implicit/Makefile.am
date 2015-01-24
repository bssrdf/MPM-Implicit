# Flag for Boost Library using Boost.m4 Macro
ACLOCAL_AMFLAGS 	= -I ./m4

AM_CPPFLAGS 		= $(BOOST_CPPFLAGS) $(CORLIB_CPPFLAGS) $(TCLAP_CPPFLAGS) $(EIGEN_CPPFLAGS)
AM_LDFLAGS 		= $(BOOST_THREADS_LDFLAGS) $(CORLIB_LDFLAGS) 

# Main Program 
bin_PROGRAMS 		= mpm

# Source Path Definition

SRC			= $(top_srcdir)/cartesian/src
MATERIAL 		= $(top_srcdir)/cartesian/material

# Files Material Libarary
MATERIALPATH 	= \
$(MATERIAL)/MaterialBaseMpm.ipp		$(MATERIAL)/MaterialFactoryMpm.hpp  	\
$(MATERIAL)/MaterialContainerMpm.hpp  	$(MATERIAL)/MaterialFactoryMpm.ipp  	\
$(MATERIAL)/MaterialBaseMpm.hpp  	$(MATERIAL)/MaterialContainerMpm.ipp 	\
$(MATERIAL)/ILE/ILE.hpp			$(MATERIAL)/ILE/ILE.ipp 		\
$(MATERIAL)/Rigid/Rigid.hpp		$(MATERIAL)/Rigid/Rigid.ipp 		\
$(MATERIAL)/Bingham/Bingham.hpp 	$(MATERIAL)/Bingham/Bingham.ipp 	\
$(MATERIAL)/Modified_Bingham/Modified_Bingham.hpp 	$(MATERIAL)/Modified_Bingham/Modified_Bingham.ipp 	

# Files SRC MPM Template
mpm_SOURCES 	= \
$(SRC)/CompareFunctor.hpp  	$(SRC)/NodeMpm.hpp		$(SRC)/ParticleSoilMpm.ipp 	\
$(SRC)/ElementMpm.hpp	   	$(SRC)/NodeMpm.ipp		$(SRC)/ParticleWaterMpm.hpp 	\
$(SRC)/ElementMpm.ipp	    	$(SRC)/ParticleCloudMpm.hpp  	$(SRC)/ParticleWaterMpm.ipp 	\
$(SRC)/Enums.hpp	    	$(SRC)/ParticleCloudMpm.ipp  	$(SRC)/SubMeshMpm.hpp		\
$(SRC)/MeshMpm.hpp	   	$(SRC)/ParticleMpm.hpp	  	$(SRC)/SubMeshMpm.ipp 		\
$(SRC)/MeshMpm.ipp	   	$(SRC)/ParticleMpm.ipp 		$(SRC)/SolverMpm.hpp		\
$(SRC)/MpmItems.hpp	   	$(SRC)/ParticleSoilMpm.hpp	$(SRC)/SolverMpm.ipp		\
$(SRC)/MatrixMpm.hpp		$(SRC)/MatrixMpm.ipp						


mpm_SOURCES 	+= 	$(MATERIALPATH) \
			$(SRC)/main.cpp 

# Define an independent executable script for inclusion in the distribution
# archive. It will not be installed on an end user's system, however.
# dist_noinst_SCRIPTS = autogen.sh
