try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

dycap_twophasenewtestBrooksCorey2l2d_pvd = PVDReader( FileName='/home/pavel/DUNE/dune-dycap/appl/twophase_test/dycap_twophase-newtest-BrooksCorey-2l2d.pvd' )

AnimationScene1 = GetAnimationScene()
AnimationScene1.EndTime = 8000.0
AnimationScene1.PlayMode = 'Snap To TimeSteps'

DataRepresentation2 = Show()

RenderView1 = GetRenderView()
a1_p_l_PVLookupTable = GetLookupTableForArray( "p_l", 1 )

DataRepresentation2.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation2.ScalarOpacityFunction = []
DataRepresentation2.ColorArrayName = 'p_l'
DataRepresentation2.ScalarOpacityUnitDistance = 0.033860043578615474
DataRepresentation2.LookupTable = a1_p_l_PVLookupTable

RenderView1.CameraViewUp = [-0.9632465915990949, 0.05823264833009091, -0.2622307427461367]
RenderView1.CameraPosition = [-0.05902948865247683, -0.4066191278790568, 1.5514334946058876]
RenderView1.CameraClippingRange = [1.3594313977853638, 2.1958953748622663]
RenderView1.CameraFocalPoint = [0.4000000059604645, 0.20000000298023224, 0.0]
RenderView1.CameraParallelScale = 0.44721360216395983
RenderView1.CenterOfRotation = [0.4000000059604645, 0.20000000298023224, 0.0]

AnimationScene1.AnimationTime = 8000.0

RenderView1.ViewTime = 8000.0
RenderView1.CameraViewUp = [0.0, 1.0, 0.0]
RenderView1.CameraPosition = [0.4000000059604645, 0.20000000298023224, 1.7279006727917343]
RenderView1.CameraFocalPoint = [0.4000000059604645, 0.20000000298023224, 0.0]
RenderView1.CameraClippingRange = [1.710621666063817, 1.7538191828836103]

a3_liquidvelocity_PVLookupTable = GetLookupTableForArray( "liquid velocity", 3 )

DataRepresentation2.ScalarOpacityFunction = []
DataRepresentation2.ColorArrayName = 'liquid velocity'
DataRepresentation2.LookupTable = a3_liquidvelocity_PVLookupTable
DataRepresentation2.ColorAttributeType = 'POINT_DATA'

Render()
