# MLIP Viewer - VR Setup Guide for Meta Quest

## Issue: VR Button Grayed Out on Meta Quest

The VR button appears grayed out because WebXR requires specific setup conditions to work properly on Meta Quest devices.

## Quick Fix Steps

### 1. Enable HTTPS (Required for WebXR)

WebXR **requires** HTTPS to function. Run these commands:

```bash
# Generate SSL certificates
npm run generate-certs

# Start server with HTTPS support
npm start
```

Then access via: `https://localhost:8443/vr.html` (instead of HTTP)

### 2. Meta Quest Browser Settings

1. **Update Browser**: Ensure you're using the latest Oculus Browser
2. **Enable WebXR**: 
   - Go to `oculus://settings` or browser settings
   - Enable "WebXR" in experimental features
   - Enable "WebGL" if not already enabled
3. **Allow Permissions**: Accept any VR/camera permission prompts

### 3. Access via HTTPS

- **Local Development**: `https://localhost:8443/vr.html`
- **Network Access**: Replace `localhost` with your computer's IP address
- **Accept Security Warning**: Click "Advanced" → "Proceed to localhost" for self-signed certificates

### 4. NEW: Debug Black Screen Issues

If you can enter VR but see only black:

1. **Use Debug Console**: Click "Show Debug" button or visit `/debug-bookmarklets.html`
2. **Check Browser Console**: Use USB debugging or on-device bookmarklets
3. **Emergency Lighting**: Use the bookmarklet to add bright lights
4. **Test Cube**: Create a red test cube to verify rendering works

## Troubleshooting

### VR Button Still Grayed Out?

1. **Check Console**: Open browser developer tools (if available) and check for errors
2. **Try Desktop First**: Test on a desktop browser with WebXR support
3. **Restart Browser**: Close and reopen the Oculus Browser
4. **Check URL**: Ensure you're using `https://` not `http://`

### Black Screen in VR Mode?

**NEW Debug Features Added:**

1. **On-Device Debug Console**: 
   - Click "Show Debug" button before entering VR
   - Press F12 or Ctrl+Shift+I to toggle debug console
   - Press Ctrl+D to run scene diagnostics

2. **Debug Bookmarklets**: 
   - Visit `https://localhost:8443/debug-bookmarklets.html`
   - Bookmark the debug tools for easy access in VR
   - Use "Emergency Lighting" if scene is too dark
   - Use "Test Cube" to verify rendering works

3. **USB Console Debugging**:
   - Enable Developer Mode on Quest
   - Connect via USB to computer
   - Open `chrome://inspect` to access full DevTools

4. **Common Black Screen Causes**:
   - **No Lighting**: Scene lights disabled or too dim
   - **No Meshes**: Molecular model failed to load
   - **Camera Position**: Camera positioned incorrectly
   - **Material Issues**: Transparent or missing materials

### Common Error Messages

- **"HTTPS required"**: You're using HTTP instead of HTTPS
- **"WebXR not supported"**: Browser doesn't support WebXR or it's disabled
- **"VR hardware not detected"**: Headset connection issues or permissions denied
- **"Scene has no lights"**: Critical lighting issue causing black screen
- **"No visible meshes"**: Molecular model not loading properly

## Alternative Solutions

### Option 1: Use ngrok for HTTPS

```bash
# Install ngrok
npm install -g ngrok

# Run your server
npm start

# In another terminal, expose with HTTPS
ngrok http 3000
```

Then use the ngrok HTTPS URL on your Quest.

### Option 2: Network HTTPS Access

1. Find your computer's IP address: `ip addr show` or `ifconfig`
2. Generate certificates for your IP:
   ```bash
   openssl req -new -x509 -keyout server-key.pem -out server-cert.pem -days 365 -nodes -subj "/CN=YOUR_IP_ADDRESS"
   ```
3. Access via `https://YOUR_IP_ADDRESS:8443/vr.html`

## What Changed

The following fixes were implemented:

1. **Enhanced HTTPS Support**: Added HTTPS server with proper headers for WebXR
2. **Better VR Detection**: Improved WebXR support detection for Meta Quest browsers
3. **Detailed Error Messages**: More informative feedback when VR fails to initialize
4. **Quest-Specific Features**: Added Meta Quest browser detection and setup instructions
5. **Fallback Handling**: Better graceful degradation when VR isn't available

## Server Features

- **Dual Mode**: Serves both HTTP (port 3000) and HTTPS (port 8443)
- **WebXR Headers**: Includes required CORS and security headers
- **Self-Signed Certs**: Automatic generation for local development
- **Quest Detection**: Recognizes Meta Quest browsers and provides specific guidance

## Production Deployment

For production use:
1. Use proper SSL certificates from a trusted CA (Let's Encrypt, etc.)
2. Set up proper domain with HTTPS
3. Ensure WebXR permissions are handled correctly
4. Test across different VR browsers and devices

## Testing VR Features

Once VR is working:
- **Controllers**: Point and click to select bonds
- **Rotation**: Hold trigger and rotate controller to modify molecular bonds
- **UI**: Use squeeze button to toggle VR interface
- **Sensitivity**: Use thumbstick to adjust rotation sensitivity

For any issues, check the browser console for detailed error messages.
