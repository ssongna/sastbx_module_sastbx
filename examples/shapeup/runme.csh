#!/bin/csh -f
                   source /home/useraccts/mcf-staff/hgliu/release/sastbx_0.99/Build/setpaths_all.csh
                   echo '<html><pre><META HTTP-EQUIV="REFRESH" CONTENT="10">' > shapeup_log.html
                   echo '<html><pre><META HTTP-EQUIV="REFRESH" CONTENT="10"> Retriving... This page will be refreshed every 10 seconds</pre><html>' > models.html
                    
sastbx.shapeup target=saxs.dat pdb=model.pdb >>& shapeup_log.html &
