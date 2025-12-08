#!/usr/bin/env python3
"""
Create a standalone HTML package from KEGG pathway markdown with all images included.
This allows the HTML to be viewed without access to the full repository.
"""

import re
import shutil
from pathlib import Path
from zipfile import ZipFile
import markdown
from markdown.extensions import codehilite, fenced_code, tables, toc

def convert_markdown_to_html(md_file, html_file):
    """Convert markdown file to HTML with proper styling and MathJax support."""
    print(f"Converting {md_file.name} to HTML...")
    
    # Read markdown
    with open(md_file, 'r', encoding='utf-8') as f:
        md_content = f.read()
    
    # Protect math equations before markdown processing
    # Replace math blocks with placeholders
    math_blocks = []
    math_counter = 0
    
    # Pattern 1: Display math ($$...$$)
    def replace_display_math(match):
        nonlocal math_counter
        math_content = match.group(1)
        placeholder = f"__MATH_DISPLAY_{math_counter}__"
        math_blocks.append(('display', math_content))
        math_counter += 1
        return placeholder
    
    # Pattern 2: Inline math ($...$)
    def replace_inline_math(match):
        nonlocal math_counter
        math_content = match.group(1)
        # Skip if it's part of a display math block
        if '$$' in math_content:
            return match.group(0)
        placeholder = f"__MATH_INLINE_{math_counter}__"
        math_blocks.append(('inline', math_content))
        math_counter += 1
        return placeholder
    
    # First, protect display math ($$...$$)
    display_math_pattern = r'\$\$([^$]+)\$\$'
    protected_content = re.sub(display_math_pattern, replace_display_math, md_content, flags=re.DOTALL)
    
    # Then, protect inline math ($...$) - but be careful not to match display math
    # Match $...$ that are not part of $$
    inline_math_pattern = r'(?<!\$)\$([^$\n]+)\$(?!\$)'
    protected_content = re.sub(inline_math_pattern, replace_inline_math, protected_content)
    
    # Configure markdown extensions
    md = markdown.Markdown(
        extensions=[
            'codehilite',
            'fenced_code',
            'tables',
            'toc',
            'attr_list',
            'def_list',
        ],
        extension_configs={
            'codehilite': {
                'css_class': 'highlight',
                'use_pygments': False,
            },
            'toc': {
                'permalink': True,
            }
        }
    )
    
    # Convert to HTML
    html_body = md.convert(protected_content)
    
    # Restore math equations
    math_counter = 0
    for math_type, math_content in math_blocks:
        if math_type == 'display':
            placeholder = f"__MATH_DISPLAY_{math_counter}__"
            # Wrap in proper div for MathJax
            math_html = f'\n<div class="math-display">\n$${math_content}$$\n</div>\n'
            html_body = html_body.replace(placeholder, math_html)
        else:  # inline
            placeholder = f"__MATH_INLINE_{math_counter}__"
            math_html = f'${math_content}$'
            html_body = html_body.replace(placeholder, math_html)
        math_counter += 1
    
    # Create full HTML document with MathJax support
    html_template = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Methods, Steps, and Results: KEGG Pathway Analysis</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            line-height: 1.6;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            color: #333;
            background-color: #fff;
        }}
        h1, h2, h3, h4, h5, h6 {{
            color: #2c3e50;
            margin-top: 1.5em;
            margin-bottom: 0.5em;
        }}
        h1 {{
            border-bottom: 3px solid #3498db;
            padding-bottom: 0.3em;
        }}
        h2 {{
            border-bottom: 2px solid #ecf0f1;
            padding-bottom: 0.3em;
        }}
        img {{
            max-width: 100%;
            height: auto;
            display: block;
            margin: 20px auto;
            border: 1px solid #ddd;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        table {{
            border-collapse: collapse;
            width: 100%;
            margin: 20px 0;
        }}
        th, td {{
            border: 1px solid #ddd;
            padding: 8px 12px;
            text-align: left;
        }}
        th {{
            background-color: #3498db;
            color: white;
            font-weight: bold;
        }}
        tr:nth-child(even) {{
            background-color: #f9f9f9;
        }}
        code {{
            background-color: #f4f4f4;
            padding: 2px 6px;
            border-radius: 3px;
            font-family: 'Courier New', monospace;
            font-size: 0.9em;
        }}
        pre {{
            background-color: #f4f4f4;
            padding: 15px;
            border-radius: 5px;
            overflow-x: auto;
            border-left: 4px solid #3498db;
        }}
        blockquote {{
            border-left: 4px solid #3498db;
            margin: 20px 0;
            padding-left: 20px;
            color: #666;
            font-style: italic;
        }}
        .math-display {{
            margin: 1.5em 0;
            text-align: center;
            overflow-x: auto;
        }}
        .toc {{
            background-color: #f9f9f9;
            border: 1px solid #ddd;
            border-radius: 5px;
            padding: 20px;
            margin: 20px 0;
        }}
        .toc ul {{
            list-style-type: none;
            padding-left: 20px;
        }}
        .toc li {{
            margin: 5px 0;
        }}
        .toc a {{
            text-decoration: none;
            color: #3498db;
        }}
        .toc a:hover {{
            text-decoration: underline;
        }}
        hr {{
            border: none;
            border-top: 2px solid #ecf0f1;
            margin: 30px 0;
        }}
    </style>
    <!-- MathJax Configuration -->
    <script>
        MathJax = {{
            tex: {{
                inlineMath: [['$', '$'], ['\\\\(', '\\\\)']],
                displayMath: [['$$', '$$'], ['\\\\[', '\\\\]']],
                processEscapes: true,
                processEnvironments: true
            }},
            options: {{
                skipHtmlTags: ['script', 'noscript', 'style', 'textarea', 'pre']
            }}
        }};
    </script>
    <!-- MathJax Library -->
    <script src="https://polyfill.io/v3/polyfill.min.js?features=es6"></script>
    <script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
</head>
<body>
{body}
</body>
</html>"""
    
    html_content = html_template.format(body=html_body)
    
    # Write HTML file
    with open(html_file, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    print(f"✓ Created HTML: {html_file}")
    return html_file

def create_standalone(md_file, html_file):
    """Create a standalone HTML package with all images."""
    base_dir = Path(__file__).parent
    project_root = base_dir.parent.parent
    standalone_dir = base_dir / 'standalone_kegg_pathway'
    images_dir = standalone_dir / 'images'
    
    # Create directories
    standalone_dir.mkdir(exist_ok=True)
    images_dir.mkdir(exist_ok=True)
    
    print(f"\nCreating standalone package in {standalone_dir}...")
    
    # Read the HTML file
    with open(html_file, 'r', encoding='utf-8') as f:
        html_content = f.read()
    
    # Find all image references (both markdown format and HTML format)
    # Pattern 1: HTML img tags: src="path"
    image_pattern1 = r'src="([^"]+\.(?:png|pdf|jpg|jpeg|gif|svg))"'
    # Pattern 2: Markdown format: ![alt](path)
    image_pattern2 = r'!\[[^\]]*\]\(([^)]+\.(?:png|pdf|jpg|jpeg|gif|svg))\)'
    
    image_matches1 = re.findall(image_pattern1, html_content, re.IGNORECASE)
    image_matches2 = re.findall(image_pattern2, html_content, re.IGNORECASE)
    all_image_paths = list(set(image_matches1 + image_matches2))
    
    print(f"Found {len(all_image_paths)} unique image references")
    
    # Track copied images to avoid duplicates
    copied_images = {}
    image_counter = 1
    
    # Copy images and update paths
    def replace_image_path(match):
        nonlocal image_counter
        full_match = match.group(0)
        img_path = match.group(1)
        
        # Skip if already processed
        if img_path in copied_images:
            new_path = copied_images[img_path]
            return full_match.replace(img_path, new_path)
        
        # Determine source path
        if img_path.startswith('results/'):
            # Path relative to project root
            source_path = project_root / img_path
        elif img_path.startswith('../../'):
            # Relative path from markdown location (../../results/...)
            source_path = base_dir / img_path
        elif img_path.startswith('../'):
            # Relative path from HTML location
            source_path = base_dir / img_path
        elif img_path.startswith('file://'):
            # Absolute file:// path
            source_path = Path(img_path.replace('file://', ''))
        else:
            # Assume relative to project root
            source_path = project_root / img_path
        
        # Prefer PNG over PDF for web display
        if source_path.suffix.lower() == '.pdf':
            png_path = source_path.with_suffix('.png')
            if png_path.exists():
                source_path = png_path
                print(f"  → Using PNG version instead of PDF: {png_path.name}")
        
        if not source_path.exists():
            print(f"  ⚠ Warning: Image not found: {source_path}")
            return full_match
        
        # Create a clean filename
        original_name = source_path.name
        # Create a unique name if needed (handle duplicates)
        if original_name in [v.split('/')[-1] for v in copied_images.values()]:
            name_parts = original_name.rsplit('.', 1)
            if len(name_parts) == 2:
                new_name = f"{name_parts[0]}_{image_counter}.{name_parts[1]}"
            else:
                new_name = f"{original_name}_{image_counter}"
        else:
            new_name = original_name
        
        # Copy image to images directory
        dest_path = images_dir / new_name
        try:
            shutil.copy2(source_path, dest_path)
            print(f"  ✓ Copied: {source_path.name} -> images/{new_name}")
        except Exception as e:
            print(f"  ✗ Error copying {source_path}: {e}")
            return full_match
        
        # Store mapping and update path
        new_relative_path = f"images/{new_name}"
        copied_images[img_path] = new_relative_path
        image_counter += 1
        
        return full_match.replace(img_path, new_relative_path)
    
    # Replace all image paths (both patterns)
    updated_html = re.sub(image_pattern1, replace_image_path, html_content, flags=re.IGNORECASE)
    updated_html = re.sub(image_pattern2, replace_image_path, updated_html, flags=re.IGNORECASE)
    
    # Remove the base tag since we're using relative paths now
    updated_html = re.sub(r'<base[^>]*>', '', updated_html)
    
    # Fix math equation rendering: Ensure math blocks are properly formatted for MathJax
    # Pattern 1: Remove <p> tags that directly wrap math blocks
    updated_html = re.sub(r'<p>\s*(\$\$.*?\$\$)\s*</p>', r'\n\1\n', updated_html, flags=re.DOTALL)
    
    # Pattern 2: Wrap standalone math blocks in divs for better MathJax processing
    def wrap_math_block(match):
        math_content = match.group(1)
        math_content = math_content.strip()
        return f'\n<div class="math-display">\n{math_content}\n</div>\n'
    
    updated_html = re.sub(r'(?<=[>\n])\s*(\$\$[^$]+\$\$)\s*(?=[<\n])', wrap_math_block, updated_html, flags=re.DOTALL)
    
    # Clean up any excessive whitespace
    updated_html = re.sub(r'\n{3,}', '\n\n', updated_html)
    
    # Write updated HTML
    standalone_html = standalone_dir / 'methods_steps_and_results_kegg_pathway.html'
    with open(standalone_html, 'w', encoding='utf-8') as f:
        f.write(updated_html)
    
    print(f"\n✓ Created standalone HTML: {standalone_html}")
    print(f"✓ Copied {len(copied_images)} unique images to {images_dir}")
    
    # Create zip file
    zip_file = base_dir / 'methods_steps_and_results_kegg_pathway_standalone.zip'
    print(f"\nCreating zip file: {zip_file}")
    
    with ZipFile(zip_file, 'w') as zipf:
        # Add HTML file
        zipf.write(standalone_html, standalone_html.name)
        print(f"  ✓ Added {standalone_html.name}")
        
        # Add all images
        for img_file in images_dir.glob('*'):
            if img_file.is_file():
                zipf.write(img_file, f"images/{img_file.name}")
                print(f"  ✓ Added images/{img_file.name}")
    
    print(f"\n✓ Standalone package created successfully!")
    print(f"  HTML: {standalone_html}")
    print(f"  Images: {images_dir} ({len(list(images_dir.glob('*')))} files)")
    print(f"  Zip: {zip_file}")
    print(f"\nTo use: Extract the zip file and open {standalone_html.name} in a browser.")

if __name__ == '__main__':
    base_dir = Path(__file__).parent
    md_file = base_dir / 'methods_steps_and_results_kegg_pathway.md'
    html_file = base_dir / 'methods_steps_and_results_kegg_pathway.html'
    
    if not md_file.exists():
        print(f"Error: Markdown file not found: {md_file}")
        exit(1)
    
    # Step 1: Convert markdown to HTML
    convert_markdown_to_html(md_file, html_file)
    
    # Step 2: Create standalone package
    create_standalone(md_file, html_file)

