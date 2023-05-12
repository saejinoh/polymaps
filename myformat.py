import streamlit as st

# ===== Settings

# ===== Utilities
def persist_widget(widget,*args,key=None,val0=None,action=lambda: None,**kwargs):
    """Persist a widget's state by using a dummy key and a callback to ensure syncing of values.

    Args:
        widget (_type_): _description_
        key (_type_): _description_
        val0 (_type_): _description_
        action (_type_): _description_

    Returns:
        _type_: _description_
    
    Notes:
        syntax is minimally invasive:
        ```
        persist_widget(widget, *(args as would be usually used), **(kwargs) as would be usually used )
        ```
        in addition, must provide the mandatory `key` and `val0` arguments explicitly by name instead of positionally.
    """
    if key is None:
        raise ValueError("must provide a key to persist the widget")
    if val0 is None:
        raise ValueError("must provide initial value of widget, via val0 argument")

    tmp_key = key + "_"
    # set up keys
    if key not in st.session_state:
        if tmp_key not in st.session_state:
            st.session_state[tmp_key] = val0
            st.session_state[key] = val0
        else:
            st.session_state[key] = st.session_state[tmp_key]
    else:
        st.session_state[tmp_key] = st.session_state[key]

    # set up callback() and return
    if widget in [st.button, st.download_button, st.form_submit_button]:
        # technically shouldn't ever be persisting a button...
        raise ValueError("Does it really make sense to persist a button value?")
        #def callback():
        #    st.session_state[key] = st.session_state[tmp_key]
        #    action()
        #return widget(*args,key=tmp_key,on_click=callback,**kwargs)
    else: #should be `on_change``
        action = kwargs["on_change"]
        def callback():
            st.session_state[key] = st.session_state[tmp_key]
            action()
        kwargs["on_change"] = callback
        return widget(*args,key=tmp_key,**kwargs)

# ===== Formatting
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

