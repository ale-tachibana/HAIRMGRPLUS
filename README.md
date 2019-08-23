# HAIRMGRPLUS
Hair Plugin for Blender - WIP

Adds a panel under particles for extra functionalities.

![HairManagerPlus](https://user-images.githubusercontent.com/30579454/63556644-07ac4f80-c51c-11e9-8bc2-6edcc88e8b12.png)

- Toggle Enable Advanced Hair

If you want to enable advanced hair AFTER you edited anything in edit mode.
  
- Copy and Paste *(still has issues)*

Makes it easier to copy and paste hair parameters between hair systems. The settings are put into the clipboard.
 
- Import Hair *(DOES NOT WORK)*

Not sure if it can really be done. The current implementation is the third attempt at doing this. The main problem is, you can't set the final calculated coordinates of each hair key. Instead you can only set the local coordinates which in theory could work if you could find the reverse transformation matrix.
  
- Export Hair *(meh works)*

The current functionality is not very useful as it can be already done by converting the hair to poly, then the poly to curve.
