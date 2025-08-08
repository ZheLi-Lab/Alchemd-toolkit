import base64
import copy
import os
import time
from .static_pics import get_cycle, get_formula
from .fast_pdf import PDFGenerator

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer, PageBreak, Image
from reportlab.lib.units import inch
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.pdfgen import canvas



cycle = get_cycle()
formula = get_formula()

collect_font = canvas.Canvas("font_check.pdf", pagesize=letter)
avail_fonts = collect_font.getAvailableFonts()
PDF_font = 'Helvetica'
if PDF_font not in avail_fonts:
    PDF_font = 'Courier'
PDF_bold_font = f"{PDF_font}-Bold"


def pdf_title_style():
    _styles = getSampleStyleSheet()
    _style = _styles['Heading2']
    _style.alignment = 0
    _style.fontSize = 20
    _style.fontName = PDF_bold_font
    _style.firstLineIndent = 0
    _style.spaceAfter = 20
    return _style


def pdf_table_style():
    table_style = [
        ('BACKGROUND', (0, 0), (-1, 0), (212 / 255, 231 / 255, 239 / 255)),
        ('TEXTCOLOR', (0, 0), (-1, 0), (0, 0, 0)),
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
        ('FONTNAME', (0, 0), (-1, 0), PDF_bold_font),
        ('FONTNAME', (0, 1), (-1, -1), PDF_font),
        ('FONTSIZE', (0, 0), (-1, -1), 12),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 5),
        ('BACKGROUND', (0, 1), (-1, -1), colors.white),
        ('GRID', (0, 1), (-1, -1), 1, colors.lightgrey),
        ('GRID', (0, 0), (-1, 0), 1, colors.black),
    ]
    return table_style


PDF_TITLE_STYLE = pdf_title_style()
TABLE_STYLE = pdf_table_style()


def load_pic(base64_pic, name, width=None, height=None):
    tmp = open(f"{name}.png", "wb+")
    tmp.write(base64.b64decode(base64_pic))
    tmp.close()
    if width is None and height is None:
        print_img = Image(f"{name}.png")
    else:
        print_img = Image(f"{name}.png", width=width, height=height, kind='proportional')
    return print_img


def static_info_page() -> list:
    _styles = getSampleStyleSheet()
    img_title_style = _styles['BodyText']
    img_title_style.fontName = PDF_font
    img_title_style.alignment = 1
    img_title_style.spaceBefore = 30
    img_title_style.spaceAfter = 30

    page_title = Paragraph("Thermodynamic Cycle of CS-FEP (Combined-Structure Free Energy Perturbation)",
                           style=PDF_TITLE_STYLE)
    cycle_pic = load_pic(cycle, name='cyc', width=400, height=400)
    formula_pic = load_pic(formula, name='formula', width=400, height=400)
    img_title = Paragraph("The thermodynamic cycle used in the CS-FEP approach.", style=img_title_style)

    return_page = [page_title, cycle_pic, img_title, formula_pic, PageBreak()]
    return return_page


def delete_pics():
    os.remove("cyc.png")
    os.remove("formula.png")


def cover_page(pair_num, ligand_num):
    _styles = getSampleStyleSheet()
    cover_title_style = _styles['Heading1']
    cover_title_style.fontName = PDF_bold_font
    cover_title_style.alignment = 1
    cover_title_style.spaceBefore = 300
    cover_title_style.spaceAfter = 100
    cover_title_style.fontSize = 60
    cover_title_style.leading = 60

    cover_desc_style = _styles['Heading2']
    cover_desc_style.alignment = 1
    cover_desc_style.fontSize = 20
    cover_desc_style.fontName = PDF_font

    sum_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())

    title_paragraph = Paragraph("Alchemd Results", style=cover_title_style)
    time_p = Paragraph(f"Date: {sum_time}", style=cover_desc_style)
    pair_p = Paragraph(f"Pairs: {pair_num}", style=cover_desc_style)
    ligand_p = Paragraph(f"Ligands: {ligand_num}", style=cover_desc_style)
    return_page = [title_paragraph, time_p, pair_p, ligand_p, PageBreak()]
    return return_page


def del_tmp_png():
    os.remove("cyc.png")
    os.remove("formula.png")


if __name__ == "__main__":
    # page = static_info_page()
    # page = cover_page(24, 16)
    # c = SimpleDocTemplate('./static_page.pdf', pagesize=letter)
    # c.build(page)
    # del_tmp_png()
    c = canvas.Canvas("table_example.pdf", pagesize=letter)
    table_data = [
        ['Edge', 'Free\nenergy\n(kcal/mol)', 'Simulation\ntime\n(ps)', 'Wall\ntime'],
        ['lig_CStoA', '0.333', '5000', '1h59m35s'],
        ['lig_CStoB', '0.333', '5000', '1h59m35s'],
        ['com_CStoA', '0.333', '5000', '1h59m35s'],
        ['com_CStoB', '0.333', '5000', '1h59m35s'],
        ['Total', '0.333', '5000', '1h59m35s'],
    ]

    PDFGenerator.draw_table(c, x=280, y=140, width=400, height=180, data=table_data, style=TABLE_STYLE)

    c.save()
