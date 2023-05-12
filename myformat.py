import streamlit as st


# Resources:
# https://discuss.streamlit.io/t/how-to-change-font-size-of-streamlit-radio-widget-title/35945/3
# https://discuss.streamlit.io/t/no-way-to-set-focus-at-top-of-screen-on-page-reload-really/15474/13
# https://discuss.streamlit.io/t/changing-the-font-size-of-widget-text-using-markdown/12933
# https://discuss.streamlit.io/t/using-custom-fonts/14005/2
# html, body, [class*="css"]  {{
#font-size: {markdown}px;
#}}
def set_font(markdown=16,widget=16):
    """
    Notes: 
        for some reason, the markdown text size option isn't working.
    """
    contents = f"""<style>
.standard-text{{
    font-size: {markdown}px;
}}
div[class*="stRadio"] > label > div[data-testid="stMarkdownContainer"] > p {{
    font-size: {widget}px;
}}
div[class*="stTextInput"] > label > div[data-testid="stMarkdownContainer"] > p {{
    font-size: {widget}px;
}}
div[class*="stSlider"] > label > div[data-testid="stMarkdownContainer"] > p {{
    font-size: {widget}px;
}}
div[class*="stTextArea"] > label > div[data-testid="stMarkdownContainer"] > p {{
    font-size: {widget}px;
}}
div[class*="row-widget stSelectbox"] > label > div[data-testid="stMarkdownContainer"] > p {{
    font-size: {widget}px;
}}
    </style>
        """
    st.markdown(contents, unsafe_allow_html=True)

def set_sidebar(max_width=650,min_width=550):
    contents = f"""
<style>
[data-testid="stSidebar"][aria-expanded="true"]{{
    max-width: {max_width}px;
    min-width: {min_width}px;
}}
    </style>
    """
    st.markdown(contents, unsafe_allow_html=True)