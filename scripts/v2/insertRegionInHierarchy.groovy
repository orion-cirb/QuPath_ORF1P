name = getSelectedObject().getName()

getSelectedObject()
runPlugin('qupath.lib.plugins.objects.DilateAnnotationPlugin', '{"radiusMicrons":-1.0,"lineCap":"ROUND","removeInterior":false,"constrainToParent":true}')
removeObject(getSelectedObject(), true)

annotation = getAnnotationObjects().find{it.getName() == name}
insertObjects(annotation)
annotation.setLocked(true)