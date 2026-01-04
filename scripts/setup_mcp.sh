#!/bin/bash
#
# Gibson Assembly Pipeline - MCP Setup Script
#
# This script configures the Google Sheets MCP server for the
# Gibson Assembly Pipeline.
#
# Usage: ./scripts/setup_mcp.sh
#
# Author: Andrew Lu (andrew.lu.chun@gmail.com)
# Affiliation: California Institute of Technology

set -e

echo "Gibson Assembly Pipeline - MCP Setup"
echo "====================================="
echo ""

# Check for required tools
if ! command -v uvx &> /dev/null; then
    echo "Error: uvx not found. Please install uv first:"
    echo "  curl -LsSf https://astral.sh/uv/install.sh | sh"
    exit 1
fi

# Check for existing configuration
if [ -f "$HOME/.config/gibson/service-account.json" ]; then
    echo "Found existing service account configuration."
    read -p "Do you want to reconfigure? (y/N): " RECONFIGURE
    if [ "$RECONFIGURE" != "y" ] && [ "$RECONFIGURE" != "Y" ]; then
        echo "Keeping existing configuration."
        exit 0
    fi
fi

# Create config directory
mkdir -p "$HOME/.config/gibson"

echo ""
echo "Step 1: Google Cloud Service Account"
echo "-------------------------------------"
echo "You need a Google Cloud service account with access to your"
echo "Google Sheets databases."
echo ""
echo "If you don't have one:"
echo "1. Go to https://console.cloud.google.com/"
echo "2. Create a project (or select existing)"
echo "3. Enable the Google Sheets API"
echo "4. Create a service account"
echo "5. Download the JSON key file"
echo "6. Share your sheets with the service account email"
echo ""

read -p "Enter path to your service account JSON file: " SA_PATH

if [ ! -f "$SA_PATH" ]; then
    echo "Error: File not found: $SA_PATH"
    exit 1
fi

# Copy to config location
cp "$SA_PATH" "$HOME/.config/gibson/service-account.json"
chmod 600 "$HOME/.config/gibson/service-account.json"

echo ""
echo "Step 2: Google Drive Folder (Optional)"
echo "---------------------------------------"
echo "If you want to access files from a specific Google Drive folder,"
echo "enter its folder ID (from the URL). Otherwise, leave blank."
echo ""

read -p "Enter Drive folder ID (or press Enter to skip): " FOLDER_ID

echo ""
echo "Step 3: Setting Environment Variables"
echo "--------------------------------------"

# Determine shell config file
if [ -n "$ZSH_VERSION" ] || [ "$SHELL" = "/bin/zsh" ]; then
    SHELL_RC="$HOME/.zshrc"
else
    SHELL_RC="$HOME/.bashrc"
fi

# Add environment variables
echo "" >> "$SHELL_RC"
echo "# Gibson Assembly Pipeline - MCP Configuration" >> "$SHELL_RC"
echo "export GIBSON_SERVICE_ACCOUNT_PATH=\"$HOME/.config/gibson/service-account.json\"" >> "$SHELL_RC"

if [ -n "$FOLDER_ID" ]; then
    echo "export GIBSON_DRIVE_FOLDER_ID=\"$FOLDER_ID\"" >> "$SHELL_RC"
fi

echo ""
echo "Added environment variables to $SHELL_RC"
echo ""

# Export for current session
export GIBSON_SERVICE_ACCOUNT_PATH="$HOME/.config/gibson/service-account.json"
if [ -n "$FOLDER_ID" ]; then
    export GIBSON_DRIVE_FOLDER_ID="$FOLDER_ID"
fi

echo "Step 4: Testing MCP Connection"
echo "------------------------------"

# Test the MCP server
echo "Testing MCP server connection..."
if uvx mcp-google-sheets@latest --help > /dev/null 2>&1; then
    echo "✓ MCP server is available"
else
    echo "⚠ Warning: Could not verify MCP server. You may need to run:"
    echo "  uvx mcp-google-sheets@latest"
fi

echo ""
echo "====================================="
echo "Setup Complete!"
echo "====================================="
echo ""
echo "Next steps:"
echo "1. Restart your terminal or run: source $SHELL_RC"
echo "2. Start Claude Code in the project directory"
echo "3. Run /sync-data to test the connection"
echo ""
echo "If you encounter issues:"
echo "- Verify your service account has access to the sheets"
echo "- Check that the Google Sheets API is enabled"
echo "- Run Claude Code with --mcp-debug for detailed logs"
echo ""
