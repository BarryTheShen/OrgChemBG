from PIL import Image, ImageDraw, ImageFont
from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas

# Define dimensions
A4_WIDTH, A4_HEIGHT = A4  # A4 dimensions in points (1 pt = 1/72 inch)
CARD_WIDTH_MM = 87
CARD_HEIGHT_MM = 57
INNER_WIDTH_MM = 79
INNER_HEIGHT_MM = 52

# Convert mm to points (1 mm = 2.83465 pts)
CARD_WIDTH_PT = CARD_WIDTH_MM * 2.83465
CARD_HEIGHT_PT = CARD_HEIGHT_MM * 2.83465
INNER_WIDTH_PT = INNER_WIDTH_MM * 2.83465
INNER_HEIGHT_PT = INNER_HEIGHT_MM * 2.83465

# Define padding and borders
BORDER_THICKNESS_PT = 2
INNER_PADDING_PT = 10 * 2.83465  # 10mm padding

# File paths
iupac_file = "iupac_names.txt"
image_folder = "output_files/"
output_pdf = "chemical_cards.pdf"

# Load IUPAC names
with open(iupac_file, "r") as file:
    iupac_names = [line.strip() for line in file.readlines()]

# Create a PDF
c = canvas.Canvas(output_pdf, pagesize=A4)
font_path = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"  # Update to the font path on your system
font_size = 12  # Default font size

# Positioning logic
x_positions = [20, 20 + CARD_WIDTH_PT + 10, 20 + 2*(CARD_WIDTH_PT + 10)]
y_positions = [A4_HEIGHT - 20 - CARD_HEIGHT_PT, A4_HEIGHT - 20 - 2*(CARD_HEIGHT_PT + 10), A4_HEIGHT - 20 - 3*(CARD_HEIGHT_PT + 20)]

for i, name in enumerate(iupac_names):
    row = i // 3
    col = i % 3
    x = x_positions[col]
    y = y_positions[row]

    # Draw the outer border
    c.setLineWidth(BORDER_THICKNESS_PT)
    c.rect(x, y, CARD_WIDTH_PT, CARD_HEIGHT_PT)

    # Draw the inner border
    inner_x = x + (CARD_WIDTH_PT - INNER_WIDTH_PT) / 2
    inner_y = y + (CARD_HEIGHT_PT - INNER_HEIGHT_PT) / 2
    c.setLineWidth(BORDER_THICKNESS_PT)
    c.rect(inner_x, inner_y, INNER_WIDTH_PT, INNER_HEIGHT_PT)

    # Draw the name
    name_font_size = min(24, font_size * INNER_WIDTH_PT / (c.stringWidth(name, font_path, font_size) + 2 * INNER_PADDING_PT))
    c.setFont(font_path, name_font_size)
    text_width = c.stringWidth(name, font_path, name_font_size)
    c.drawString(inner_x + (INNER_WIDTH_PT - text_width) / 2, y + CARD_HEIGHT_PT - INNER_PADDING_PT - name_font_size, name)

    # Draw the image
    img = Image.open(f"{image_folder}/compound_{i+1}.png").resize((300, 300))
    img_x = inner_x + (INNER_WIDTH_PT - 300) / 2
    img_y = inner_y + (INNER_HEIGHT_PT - 300) / 2
    img.save(f"temp_img_{i}.png")
    c.drawImage(f"temp_img_{i}.png", img_x, img_y, width=300, height=300)

    # Move to the next card
    if (i + 1) % 9 == 0:
        c.showPage()

# Save the PDF
c.save()

print(f"PDF saved as {output_pdf}")
