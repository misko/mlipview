# Meta Quest Browser Debug Console Guide

## Method 1: USB Debugging (Recommended)

1. **Enable Developer Mode on Quest**:
   - Open Meta Quest app on your phone
   - Go to Settings → Your headset → Developer Mode
   - Toggle on Developer Mode
   - https://developers.meta.com/horizon/documentation/native/android/mobile-device-setup/

2. **Connect via USB**:
   - Connect Quest to computer via USB-C cable
   - Put on headset and allow USB debugging when prompted

3. **Use Chrome DevTools**:
   - Open Chrome on your computer
   - Go to `chrome://inspect`
   - Your Quest browser should appear under "Remote Target"
   - Click "Inspect" to open DevTools

## Method 2: Remote Debugging via IP

1. **Find Quest's IP Address**:
   - In Quest browser, go to `about:version`
   - Note the IP address

2. **Enable ADB on Computer**:
   ```bash
   adb connect YOUR_QUEST_IP:5555
   ```

3. **Port Forward**:
   ```bash
   adb forward tcp:9222 localabstract:chrome_devtools_remote
   ```

4. **Access DevTools**:
   - Open `localhost:9222` in your computer's browser

## Method 3: On-Device Console (Limited)

1. **URL Bar Method**:
   - In Quest browser, type: `javascript:console.log("test")`
   - Press Enter to execute

2. **Bookmarklet Console**:
   - Bookmark this JavaScript: 
   ```javascript
   javascript:(function(){var script=document.createElement('script');script.src='https://cdn.jsdelivr.net/npm/eruda';document.body.appendChild(script);script.onload=function(){eruda.init();}})();
   ```
   - Click bookmark to show on-device console

## Quick Debug Commands

Once you have console access, try these:

```javascript
// Check VR status
console.log('VR Scene:', window.scene);
console.log('VR Camera:', window.scene?.activeCamera);
console.log('VR Lights:', window.scene?.lights);

// Check if objects are loaded
console.log('Meshes:', window.scene?.meshes?.length);
console.log('Materials:', window.scene?.materials?.length);
```
