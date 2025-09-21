#!/bin/bash

# Generate self-signed SSL certificates for local HTTPS testing
# Required for WebXR functionality

echo "Generating self-signed SSL certificates for localhost..."

# Generate private key
openssl genrsa -out localhost-key.pem 2048

# Generate certificate signing request
openssl req -new -key localhost-key.pem -out localhost.csr -subj "/C=US/ST=State/L=City/O=Organization/CN=localhost"

# Generate self-signed certificate
openssl x509 -req -in localhost.csr -signkey localhost-key.pem -out localhost-cert.pem -days 365

# Clean up CSR file
rm localhost.csr

echo "SSL certificates generated:"
echo "  - localhost-key.pem (private key)"
echo "  - localhost-cert.pem (certificate)"
echo ""
echo "To use HTTPS with WebXR:"
echo "1. Run: npm start"
echo "2. Open: https://localhost:8443/vr.html"
echo "3. Accept the security warning (self-signed certificate)"
echo "4. On Meta Quest, you may need to enable 'Allow insecure content' in browser settings"
echo ""
echo "Note: For production, use proper SSL certificates from a trusted CA."
