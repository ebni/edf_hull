#include      "qhull-src/src/libqhull_r/libqhull_r.h"
#include      <limits.h>

qhT qh_qh = {
	.ALLpoints = 0,
	.ANGLEmerge = 1,             /* OK */
	.APPROXhull = 0,
	.MINoutside = 0,
	.ANNOTATEoutput = 0,
	.ATinfinity = 0,
	.AVOIDold = 0,
	.BESToutside = 0,
	.CDDinput = 0,
	.CDDoutput = 0,
	.CHECKfrequently = 0,
	.premerge_cos = REALmax,     /* OK */
	.postmerge_cos = REALmax,    /* OK */
	.DELAUNAY = 0,
	.DOintersections = 0,
	.DROPdim = -1,               /* OK */
	.FLUSHprint = 0,
	.FORCEoutput = 0,
	.GOODpoint = 0,
	.GOODpointp = 0x0,
	.GOODthreshold = 0,
	.GOODvertex = 0,
	.GOODvertexp = 0x0,
	.HALFspace = 0,
	.ISqhullQh = 0,
	.IStracing = 0,
	.KEEParea = 0,
	.KEEPcoplanar = 0,
	.KEEPinside = 0,
	.KEEPmerge = 0,
	.KEEPminArea = REALmax,      /* OK */
	.MAXcoplanar = REALmax,      /* OK */
	.MERGEexact = 0,
	.MERGEindependent = 1,       /* OK */
	.MERGING = 0,
	.premerge_centrum = 0,
	.postmerge_centrum = 0,
	.MERGEpinched = 1,           /* OK */
	.MERGEvertices = 1,          /* OK */
	.MINvisible = REALmax,       /* OK */
	.NOnarrow = 0,
	.NOnearinside = 0,
	.NOpremerge = 0,
	.ONLYgood = 0,
	.ONLYmax = 0,
	.PICKfurthest = 0,
	.POSTmerge = 0,
	.PREmerge = 0,
	.PRINTcentrums = 0,
	.PRINTcoplanar = 0,
	.PRINTdim = 0,
	.PRINTdots = 0,
	.PRINTgood = 0,
	.PRINTinner = 0,
	.PRINTneighbors = 0,
	.PRINTnoplanes = 0,
	.PRINToptions1st = 0,
	.PRINTouter = 0,
	.PRINTprecision = 1,         /* OK */
	.PRINTout = {qh_PRINTextremes, 0},    /* OK */
	.PRINTridges = 0,
	.PRINTspheres = 0,
	.PRINTstatistics = 0,
	.PRINTsummary = 0,
	.PRINTtransparent = 0,
	.PROJECTdelaunay = 0,
	.PROJECTinput = 0,
	.RANDOMdist = 0,
	.RANDOMfactor = 0,
	.RANDOMa = 0,
	.RANDOMb = 0,
	.RANDOMoutside = 0,
	.REPORTfreq = 0,
	.REPORTfreq2 = 0,
	.RERUN = 0,
	.ROTATErandom = INT_MIN,     /* OK */
	.SCALEinput = 0,
	.SCALElast = 0,
	.SETroundoff = 0,
	.SKIPcheckmax = 0,
	.SKIPconvex = 0,
	.SPLITthresholds = 0,
	.STOPcone = 0,
	.STOPpoint = 0,
	.TESTpoints = 0,
	.TESTvneighbors = 0,
	.TRACElevel = 0,
	.TRACElastrun = 0,
	.TRACEpoint = qh_IDunknown,    /* OK */
	.TRACEdist = REALmax,          /* OK */
	.TRACEmerge = 0,
	.TRIangulate = 0,
	.TRInormals = 0,
	.UPPERdelaunay = 0,
	.USEstdout = 0,
	.VERIFYoutput = 0,
	.VIRTUALmemory = 0,
	.VORONOI = 0,
	.AREAfactor = 0,
	.DOcheckmax = 0,
	.feasible_string = 0x0,
	.feasible_point = 0x0,
	.GETarea = 0,
	.KEEPnearinside = 0,
	.hull_dim = 0,      /* src/libqhull_r/global_r.c:1507 */  /* TO BE SET TO THE DIMENSION OF POINTS: 3 in the example */
	.input_dim = 0,     /* src/libqhull_r/global_r.c:1507 */  /* TO BE SET TO THE DIMENSION OF POINTS: 3 in the example */
	.num_points = 0,    /* src/libqhull_r/global_r.c:1506 */  /* TO BE SET TO THE NUMBER OF POINTS: 126 in the example */
	.first_point = 0x0, /* src/libqhull_r/global_r.c:1505 */  /* TO BE SET TO THE ADDRESS OF THE POINTS */
	.POINTSmalloc = 0, /* src/libqhull_r/global_r.c:1504 */
	.input_points = 0x0,
	.input_malloc = 0,
	.qhull_command = "qhull Fx TI data_in.txt TO data_out.txt",   /* OK */   /* FORSE SI PUO` ANCHE TOGLIERE */
	.qhull_commandsiz2 = 0,
	.rbox_command = {0},
	.qhull_options = "  run-id 772309079  Fxtremes  TInput-file  data_in.txt  TOutput-file\n  data_out.txt",
	.qhull_optionlen = 14,      /* SE IGNORABILE, MEGLIO*/
	.qhull_optionsiz = 0,
	.qhull_optionsiz2 = 0,
	.VERTEXneighbors = 0,
	.ZEROcentrum = 0,
	.upper_threshold = 0x0,
	.lower_threshold = 0x0,
	.upper_bound = 0x0,
	.lower_bound = 0x0,
	.ANGLEround = 0,
	.centrum_radius = 0,
	.cos_max = 0,
	.DISTround = 0,
	.MAXabs_coord = 0,
	.MAXlastcoord = 0,
	.MAXsumcoord = 0,
	.MAXwidth = -REALmax,   /* OK */
	.MINdenom = 0,
	.MINdenom_1_2 = 0,
	.MINdenom_2 = 0,
	.MINlastcoord = 0,
	.NARROWhull = 0,
	.NEARzero = 0x0,
	.NEARinside = 0,
	.ONEmerge = 0,
	.outside_err = REALmax,  /* OK */
	.WIDEfacet = 0,
	.qhull = "qhull",
	.jmpXtra = {0},
	.restartexit = {{.__jmpbuf = {0}, .__mask_was_saved = 0, .__saved_mask = {.__val = {0}}}},
	.jmpXtra2 = {0},
	.interior_point = 0x0,
	.center_size = 0,
	.TEMPsize = 0,
	.facet_list = 0x0,
	.facet_tail = 0x0,
	.facet_next = 0x0,
	.newfacet_list = 0x0,
	.visible_list = 0x0,
	.num_visible = 0,
	.tracefacet_id = UINT_MAX,   /* OK */
	.tracefacet = 0x0,
	.tracevertex_id = UINT_MAX,  /* OK */
	.tracevertex = 0x0,
	.vertex_list = 0x0,
	.vertex_tail = 0x0,
	.newvertex_list = 0x0,
	.num_facets = 0,
	.num_vertices = 0,
	.num_outside = 0,
	.num_good = 0,
	.facet_id = 0,
	.ridge_id = 0,
	.vertex_id = 0,
	.first_newfacet = 0,
	.hulltime = 0,
	.ALLOWrestart = 0,
	.build_cnt = 0,
	.CENTERtype = qh_ASnone,
	.furthest_id = qh_IDunknown,  /* OK */
	.GOODclosest = 0x0,
	.hasAreaVolume = 0,
	.hasTriangulation = 0,
	.isRenameVertex = 0,
	.JOGGLEmax = REALmax,  /* OK */
	.maxoutdone = 0,
	.max_outside = 0,
	.max_vertex = 0,
	.min_vertex = 0,
	.NEWfacets = 0,
	.NEWtentative = 0,
	.findbestnew = 0,
	.findbest_notsharp = 0,
	.NOerrexit = 0,
	.PRINTcradius = 0,
	.PRINTradius = 0,
	.POSTmerging = 0,
	.printoutvar = 0,
	.printoutnum = 0,
	.repart_facetid = 0,
	.retry_addpoint = 0,
	.QHULLfinished = 0,
	.totarea = 0,
	.totvol = 0,
	.visit_id = 0,
	.vertex_visit = 0,
	.WAScoplanar = 0,
	.ZEROall_ok = 0,
	.facet_mergeset = 0x0,
	.degen_mergeset = 0x0,
	.vertex_mergeset = 0x0,
	.hash_table = 0x0,
	.other_points = 0x0,
	.del_vertices = 0x0,
	.gm_matrix = 0x0,
	.gm_row = 0x0,
	.line = 0x0,
	.half_space = 0x0,
	.temp_malloc = 0x0,
	.ERREXITcalled = 0,
	.firstcentrum = 0,
	.old_randomdist = 0,
	.coplanarfacetset = 0x0,
	.last_low = REALmax,      /* OK */
	.last_high = REALmax,     /* OK */
	.last_newhigh = REALmax,  /* OK */
	.lastreport = 0,
	.mergereport = 0,
	.old_tempstack = 0x0,
	.ridgeoutnum = 0,
	.rbox_errexit = {{.__jmpbuf = {0, 0, 0, 0, 0, 0, 0, 0},
			  .__mask_was_saved = 0,
			  .__saved_mask = {.__val = {0}}}},
	.jmpXtra3 = {0},
	.rbox_isinteger = 0,
	.rbox_out_offset = 0,
	.cpp_object = 0x0,
	.qhmem = {.BUFsize = 0,
		  .BUFinit = 0,
		  .TABLEsize = 0,
		  .NUMsizes = 0,
		  .LASTsize = 0,
		  .ALIGNmask = 0,
		  .freelists = 0x0,
		  .sizetable = 0x0,
		  .indextable = 0x0,
		  .curbuffer = 0x0,
		  .freemem = 0x0,
		  .freesize = 0,
		  .tempstack = 0x0,
		  .IStracing = 0,
		  .cntquick = 0,
		  .cntshort = 0,
		  .cntlong = 0,
		  .freeshort = 0,
		  .freelong = 0,
		  .totbuffer = 0,
		  .totdropped = 0,
		  .totfree = 0,
		  .totlong = 0,
		  .maxlong = 0,
		  .totshort = 0,
		  .totunused = 0,
		  .cntlarger = 0,
		  .totlarger = 0
	},
};



