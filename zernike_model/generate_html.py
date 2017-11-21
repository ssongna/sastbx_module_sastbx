def generate_jmol_html(ave_maps, ave_ccs, ave_cutoffs, model_list, ccs, cutoffs, cids, html_file, pdb=None):
  from libtbx.test_utils import run_command
  plot_template='''
<head>

<script type="text/javascript" src="/openflashchart/js/swfobject.js"></script>
<script type="text/javascript">
swfobject.embedSWF("/openflashchart/open-flash-chart.swf", "my_chart", "440", "330", "9.0.0","expressInstall.swf", {"data-file":"%s"});
</script>

</head>
<body>
<div id="my_chart"></div>
'''

  template_w_pdb= '''
<tr>
 <td>
  <script>
   script0 =
   'load "%s"; cartoon only; color structure; isosurface cutoff %f "%s";background white;color isosurface translucent; center $isosurface1;zoom 100; set spiny 10;set spinfps 15;spin off;';
   jmolApplet(480, script0);
  </script>
 </td>
 <td>
 <script type="text/javascript">
   jmolButton("spacefill 1.5; color yellowtint;", "spacefill")
   jmolBr()
   jmolButton("cartoon only; color structure;", "cartoon")
 </script>
 </td>
 <td>
  <script type="text/javascript">
   swfobject.embedSWF("/openflashchart/open-flash-chart.swf", "my_chart", "480", "360", "9.0.0","expressInstall.swf", {"data-file":"%s"});
  </script>
  <div id="my_chart"></div>
 </td>
 <td>%s</td>
 <td><center>%3.2f</center></td>
 <td><center>%d</center></td>
 <td><a href="http://www.rcsb.org/pdb/explore/explore.do?structureId=%s" target="_blank"> %s  </a></td>
</tr>
'''  # Need pdb_name, cutoff, map_name, map_id, cc, cluster_id, pdb_id, pdb_id

  template_wo_pdb= '''
<tr>
 <td>
  <script>
   script0 =
   'isosurface cutoff %f "%s";background white;color isosurface translucent; center $isosurface1;zoom 30; set spiny 10;set spinfps 15;spin off;';
   jmolApplet(480, script0);
  </script>
 </td>
 <td>
  <script type="text/javascript">
   swfobject.embedSWF("/openflashchart/open-flash-chart.swf", "my_chart", "480", "360", "9.0.0","expressInstall.swf", {"data-file":"%s"});
  </script>
  <div id="my_chart"></div>
 </td>
 <td>%s</td>
 <td><center>%3.2f</center></td>
 <td><center>%d</center></td>
 <td><a href="http://www.rcsb.org/pdb/explore/explore.do?structureId=%s" target="_blank"> %s  </a></td>
  </tr>
'''  # Need cutoff, map_name, json.iq map_id, cc, cluster_id, pdb_id, pdb_id


  light_template_w_pdb= '''
<tr>
 <td>
  <script>
   script0 =
   'load "%s"; cartoon only; color structure; isosurface cutoff %f "%s";background white;color isosurface translucent; center $isosurface1;zoom 100; set spiny 10;set spinfps 15;spin off;';
   jmolApplet(300, script0);
  </script>
 </td>
 <td>
 <script type="text/javascript">
   jmolButton("spacefill 1.5; color yellowtint;", "spacefill")
   jmolBr()
   jmolButton("cartoon only; color structure;", "cartoon")
 </script>
 </td>
 <td>%s</td>
 <td><center>%3.2f</center></td>
 <td><center>%d</center></td>
 <td><center><img src="cluster.png" alt="image not generated, check neato" width="250" height="250" /></center></td>
</tr>
'''  # Need pdb_name, cutoff, map_name, map_id, cc, cluster_id

  light_template_wo_pdb= '''
<tr>
 <td>
  <script>
   script0 =
   'isosurface cutoff %f "%s";background white;color isosurface translucent; center $isosurface1;zoom 30; set spiny 10;set spinfps 15;spin off;';
   jmolApplet(300, script0);
  </script>
 </td>
 <td>%s</td>
 <td><center>%3.2f</center></td>
 <td><center>%d</center></td>
 <td><center><img src="cluster.png" alt="image not generated, check neato" width="250" height="250" /></center></td>
</tr>
'''  # Need cutoff, map_name, map_id, cc, cluster_id


  link_template='''
 <tr>
 <td><center>%d</center></td>
 <td>%s</td>
 <td><center>%3.2f</center></td>
 <td><a href="%s"> Larger View </a></td>
 <td><center>%d</center></td>
 <td><a href="http://www.rcsb.org/pdb/explore/explore.do?structureId=%s" target="_blank"> %s  </a></td>
 <td><a href=%s> %s </a> </td>
  </tr>
''' # model_name, cc, model_html, cluster_id, pdb_id, pdb_id, ccp4_file, ccp4_file

  html_head='''<html>
<head>
 <title> Models from Shapeup</title>
 <script src="/jmol/Jmol.js" type="text/javascript"></script>
</head>
[<a href="/cgi-bin/index.html">Index of Services</a>]
[<a href="/cgi-bin/shapeup.html">New input</a>]
<hr>
<body>
  <script type="text/javascript">
  jmolInitialize("/jmol");
  </script>
 <h3> Top 10 Models Found </h3>
 <hr>
 <table border="1">
 <tr>
 <th> Model </th>
 <th> Name  </th>
 <th> Correlation Coefficient<small><a href="#cc">*</a></small> </th>
 <th> Larger View </th>
 <th> cluster ID </th>
 <th> Link to PDB </th>
 <th> Link to Map (ccp4)</th>
 </tr>
''' ## Head of Html and setup table headers

  individual_html_head='''<html>
<head>
 <title> Models from Shapeup</title>
 <script src="/jmol/Jmol.js" type="text/javascript"></script>
 <script type="text/javascript" src="/openflashchart/js/swfobject.js"></script>
</head>
[<a href="/cgi-bin/index.html">Index of Services</a>]
[<a href="/cgi-bin/shapeup.html">New input</a>]
[<a href="models.html">Model Summary Page</a>]
<hr>
<body>
  <script type="text/javascript">
  jmolInitialize("/jmol");
  </script>
 <table border="1">
 <tr>
 <th> Model </th>
 <th> Representation </th>
 <th> Data  </th>
 <th> Name  </th>
 <th> Correlation Coefficient<small><a href="#cc">*</a></small> </th>
 <th> cluster ID </th>
 <th> Link to PDB</th>
 </tr>
''' ## Head of Html and setup table headers

  individual_html_head_wo_pdb='''<html>
<head>
 <title> Models from Shapeup</title>
 <script src="/jmol/Jmol.js" type="text/javascript"></script>
 <script type="text/javascript" src="/openflashchart/js/swfobject.js"></script>
</head>
[<a href="/cgi-bin/index.html">Index of Services</a>]
[<a href="/cgi-bin/shapeup.html">New input</a>]
[<a href="models.html">Model Summary Page</a>]
<hr>
<body>
  <script type="text/javascript">
  jmolInitialize("/jmol");
  </script>
 <table border="1">
 <tr>
 <th> Model </th>
 <th> Data  </th>
 <th> Name  </th>
 <th> Correlation Coefficient<small><a href="#cc">*</a></small> </th>
 <th> cluster ID </th>
 <th> Link to PDB</th>
 </tr>
''' ## Head of Html and setup table headers


  if( pdb is not None ): template = template_w_pdb
  else: template=template_wo_pdb

  out = open( html_file, 'w')
  print>>out, html_head
  ii = 1
  for cut, model, cc, cid in zip( cutoffs, model_list, ccs, cids ):
    model_id = model.split('.')[0]
    pdb_id = model_id.split('_')[1]
    model_html = pdb_id+'.html'
    json_name="query_"+str(ii)+".iq.json"
    print>>out, link_template%(ii, pdb_id, cc, model_html,cid, pdb_id, pdb_id, model, model)
    ii = ii + 1

    this_out=open(model_html, 'w')
    if( pdb is not None):
      print>>this_out, individual_html_head
      print>>this_out, template%(pdb,cut,model,json_name, pdb_id, cc, cid, pdb_id, pdb_id)
    else:
      print>>this_out, individual_html_head_wo_pdb
      print>>this_out, template%(cut,model,json_name, pdb_id, cc, cid, pdb_id, pdb_id)
    print>>this_out, '''</table>'''
    print>>this_out, '''<p><a name="cc">* Correlation Coeffient is computed w.r.t. the given PDB model (if any) or the Best Model found </a>'''
    print>>this_out, '''</body></html>'''
    this_out.close()

  print>>out, '''
 </table>
 <hr>
 <h3> Average Maps: </h3>
 <table border="1">
 <tr>
 <th> Model </th>
 <th> Representation </th>
 <th> Name  </th>
 <th> Correlation Coefficient<small><a href="#cc">*</a></small> </th>
 <th> cluster ID </th>
 <th> cluster connection </th>
 </tr>
'''

  ii = 1
  for cut, map, cc in zip( ave_cutoffs, ave_maps, ave_ccs):
    model_id="average_model_"+str(ii)
    if( pdb is not None):
      print>>out, light_template_w_pdb%(pdb,cut,map,model_id, cc, ii)
    else:
      print>>out, light_template_wo_pdb%(cut,map,model_id, cc, ii)
    ii = ii + 1

  print>>out, '''</table><hr>'''
#  run_command('/sw/bin/tar -czf models.tgz *')
  print>>out, '''Download results in tar file <a href="models.tgz"> Models</a>.'''
  print>>out, '''<p><a name="cc">* Correlation Coeffient is computed w.r.t. the given PDB model (if any) or the Best Model found </a>'''
  print>>out, '''<br>No Cartoon Model? This may help. <a href="/cgi-bin/faq.html#5"> FAQ 5</a><br>'''
  print>>out, '''</body></html>'''
  out.close()
  return

if __name__ == "__main__":
  # map_list cc_list cluster_id_list
  #def generate_jmol_html(map_files, ccs, cutoffs, cids, html_file, pdb=None):
  generate_jmol_html(["test.ccp4"],[1.0],[0.07],["test.ccp4"], [1.0], [0.07], [1], "models.html", pdb="test.pdb")
