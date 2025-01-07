# conf.py

import os
import sys
sys.path.insert(0, os.path.abspath('.'))

# -- 项目信息 -----------------------------------------------------
project = 'VeloGIF'
author = '49'
copyright = '2024, 49'

# -- 常规配置 ---------------------------------------------------
extensions = [
    # 在此添加其他需要的扩展，例如：
    # "sphinx.ext.autodoc",
    # "sphinx.ext.napoleon",
    'sphinx_copybutton',
]

templates_path = ['_templates']
exclude_patterns = []

# -- HTML 输出配置 -------------------------------------------------
html_theme = "sphinx_rtd_theme"

# 移除不支持的主题选项 'github_url'
html_theme_options = {
   
}

# 隐藏“View page source”链接
html_show_sourcelink = False

html_static_path = ['_static']
html_css_files = [
    'css/custom.css',
]
html_js_files = [
    'js/custom.js',
]
# 使用 html_context 配置 GitHub 集成
html_context = {
    'display_github': True,
    'github_user': 'liyarubio',
    'github_repo': 'VeloGIF',
    'github_version': 'main',  
    'conf_py_path': '/sphinx/source/' 
}
