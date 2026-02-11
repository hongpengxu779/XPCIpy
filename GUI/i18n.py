# -*- coding: utf-8 -*-
"""
Centralized Chinese translations for the XPCIpy GUI.
All UI-facing strings are stored here for easy maintenance.
"""

# ==================== Main Window ====================
STATUS_READY = "就绪"
MENU_HELP = "帮助"
MENU_HELP_GUIDE = "帮助 / 用户指南"
MENU_CITE = "如何引用 XPCIpy"
MENU_LICENSE = "许可证"

# ==================== Loading Screen ====================
LOADING_TEXT = "正在加载 XPCIpy…"

# ==================== Tab Names ====================
TAB_TLREC = "TL重建"
TAB_INLINE = "在线仿真"
TAB_CHECK_TL = "检查Talbot-Lau效应"
TAB_TL_SIM = "Talbot-Lau相衬仿真"
TAB_TLREC_BATCH = "TL批量重建(开发中)"

# ==================== Common Labels ====================
LBL_N_PIXELS = "n: 波前像素大小"
LBL_PIXEL_SIZE = "像素尺寸 (微米)"
LBL_PIXEL_SIZE_MICRON = "像素尺寸 (微米)"
LBL_SOURCE_FWHM = "源 FWHM (微米)"
LBL_BEAM_SHAPE = "光束形状"
LBL_SPECTRUM = "能谱"
LBL_ENERGY = "能量 (keV)"
LBL_DESIGN_ENERGY = "设计能量 (keV)"
LBL_OBJECT = "物体"
LBL_IMAGE = "图像"
LBL_DETECTOR_PX = "探测器像素尺寸 (微米)"
LBL_DETECTOR_RES_FWHM = "探测器分辨率 (FWHM 微米)"
LBL_DETECTOR_RES_PX = "探测器分辨率 (像素尺寸 微米)"

# ==================== Inline Simulation Labels ====================
LBL_DOD = "物体到探测器距离 (cm)"
LBL_DSO = "源到物体距离 (cm)"

# ==================== Check TL Labels ====================
LBL_GRATING_PERIOD = "光栅周期 (微米)"
LBL_DUTY_CYCLE = "占空比"
LBL_MATERIAL_CUSTOM = "材料 (仅用于自定义光栅)"
LBL_BAR_HEIGHT = "栅条高度 (微米, 仅自定义光栅)"
LBL_GRATING_TYPE = "光栅类型"
LBL_TALBOT_MULTIPLES = "Talbot距离倍数"
LBL_CALC_COUNT = "计算次数"
LBL_TALBOT_DIST = "Talbot距离 (cm)"

# ==================== TL Simulation Labels ====================
LBL_DSG1 = "源到G1距离 (cm)"
LBL_DOG1 = "物体到G1距离 (cm)"
LBL_MAGNIFICATION = "放大倍数"
LBL_DG1G2 = "G1-G2距离 (cm)"
LBL_G1_PERIOD = "G1周期 (微米)"
LBL_G2_PERIOD = "G2周期 (微米)"
LBL_G1_PHASE = "G1相位"
LBL_MOVABLE_GRATING = "可移动光栅"
LBL_NUM_STEPS = "步数 (整数)"
LBL_STEP_LENGTH = "步长 (微米)"

# ==================== Buttons ====================
BTN_RUN = "运行"
BTN_EXIT = "退出"
BTN_SAVE_PRESET = "保存预设"
BTN_LOAD_PRESET = "加载预设"
BTN_SET_OBJ_PARAMS = "设置物体参数"
BTN_SAVE_IMAGE = "保存图像"
BTN_SAVE_STACK_OBJ = "保存物体图像栈"
BTN_SAVE_STACK_REF = "保存参考图像栈"
BTN_SEND_TLREC = "发送到TL重建"
BTN_CREATE_ZIP = "创建仿真数据ZIP文件"
BTN_APPLY_CHANGES = "应用更改"
BTN_SAVE_MODIFY = "保存/修改"
BTN_MODIFY_PARAMS = "修改默认参数"
BTN_RETRIEVE = "重建"
BTN_UPDATE_MC = "更新调制曲线"
BTN_APPLY = "应用"

# ==================== Object Parameters Window ====================
WIN_OBJ_PARAMS = "物体参数"
LBL_RADIUS = "半径 (微米):"
LBL_OUTER_RADIUS = "外半径 (微米):"
LBL_INNER_RADIUS = "内半径 (微米):"
LBL_XSHIFT = "X方向偏移 (像素):"
LBL_YSHIFT = "Y方向偏移 (像素):"
LBL_ORIENTATION = "方向"
LBL_MATERIAL = "材料"

# ==================== TLRec GUI ====================
WIN_DEFAULT_PARAMS = "默认参数"
LBL_G1_PERIOD_UM = "G1周期 (微米):"
LBL_G2_PERIOD_UM = "G2周期 (微米):"
LBL_DESIGN_ENERGY_KEV = "设计能量 (keV):"
LBL_DSG1_CM = "源到G1距离 (cm):"
LBL_DOG1_CM = "物体到G1距离 (cm):"
LBL_DG1G2_CM = "G1到G2距离 (cm):"
LBL_DET_PX_UM = "探测器像素尺寸 (微米):"
LBL_FREQ_CUTOFF = "频率截止 (维纳滤波器):"
LBL_BUTTERWORTH_N = "巴特沃斯滤波器阶数 (整数):"
LBL_BUTTERWORTH_S = "巴特沃斯滤波器信号:"
LBL_ALGORITHM = "算法:"
LBL_REC_ALGORITHM = "重建算法"
LBL_X_COORD = "x 坐标"
LBL_Y_COORD = "y 坐标"
FRAME_RECONSTRUCTION = "重建"

# ==================== DropZone ====================
DZ_REF_TITLE = "参考图像栈"
DZ_REF_DESC = "将TIFF图像栈拖放到此处"
DZ_REF_BTN = "浏览参考图像"
DZ_OBJ_TITLE = "物体图像栈"
DZ_OBJ_DESC = "将TIFF图像栈拖放到此处"
DZ_OBJ_BTN = "浏览物体图像"
DZ_NO_FILE = "未加载文件"
DZ_LOADED = "已加载: {}"

# ==================== Status Messages ====================
ST_RUNNING_INLINE = "正在运行在线仿真..."
ST_RUNNING_TL = "正在运行Talbot-Lau仿真..."
ST_RUNNING_CARPET = "正在运行Talbot毯仿真..."
ST_ERR_INLINE = "在线仿真的物理参数有误。"
ST_ERR_CARPET = "Talbot毯仿真的物理参数有误。"
ST_ERR_TL = "Talbot-Lau仿真的物理参数有误。"
ST_PRESET_SAVED = "预设保存成功。"
ST_PRESET_LOADED = "预设加载成功。"
ST_REF_LOADED = "参考图像已加载。请加载物体图像。"
ST_OBJ_LOADED = "物体图像已加载，请加载参考图像。"
ST_IMAGES_LOADED = "图像加载完成！"
ST_LOADING_REF = "正在加载参考图像..."
ST_LOADING_OBJ = "正在加载物体图像..."
ST_PREPARING_REC = "正在准备使用 {} 进行重建..."
ST_RUNNING_REC = "正在运行重建 ({})..."
ST_REC_FINISHED = "重建完成，耗时 {:.2f} 秒。"
ST_SELECTED_ALGO = "已选择重建算法: {}"
ST_SIM_LOADED = "仿真图像已加载到TL重建。"
ST_TLREC_NA = "TL重建界面不可用。"
ST_ERROR = "错误。"
ST_WIENER_APPLIED = "维纳滤波已应用 (v0={}, n={}, s={})"
ST_TL_INVALID = "Talbot配置在这些参数下物理无效。"
ST_REC_ERROR = "重建过程中出错。"

# ==================== Plot Titles ====================
PLOT_INLINE = "在线仿真"
PLOT_PHASE_STEPPING = "相位步进曲线"
PLOT_ONE_PROJ = "单次投影"
PLOT_PHASE_GRAD = "相位梯度"
PLOT_INT_PHASE = "积分相位"
PLOT_TRANSMISSION = "透射"
PLOT_DARK_FIELD = "暗场"

# ==================== Wiener Filter Panel ====================
FRAME_WIENER = "维纳滤波器"
FRAME_WL_CONTROLS = "窗宽/窗位控制"

# ==================== Validation Messages ====================
VAL_INVALID_PARAMS = "参数无效"
VAL_POSITIVE_PERIOD = "请确保G1周期、能量和Talbot倍数为正值。"
VAL_POSITIVE_FWHM = "请确保源FWHM为正值。"
VAL_SPHERE_RADIUS = "球体半径必须 > 0。"
VAL_CYL_OUTER = "圆柱体外半径必须 > 0。"
VAL_CYL_INNER_NEG = "圆柱体内半径不能为负。"
VAL_CYL_INNER_GE = "圆柱体内半径必须小于外半径。"
VAL_POSITIVE_N = "请确保像素数为正值。"
VAL_POSITIVE_DET = "请确保探测器像素尺寸和探测器FWHM为正值。"
VAL_POSITIVE_DSO = "请确保DSO和DOG1为正值。"
VAL_POSITIVE_PX = "请确保像素尺寸为正值。"
VAL_POSITIVE_DSO_DOD = "请确保DSO和DOD为正值。"
VAL_POSITIVE_ENERGY = "请确保能量为正值。"
VAL_POSITIVE_PERIOD_G = "请确保光栅周期为正值。"
VAL_DC_RANGE = "请确保占空比在0到1之间。"
VAL_CUSTOM_BAR = "自定义光栅的栅条高度必须 > 0 微米。"
VAL_POSITIVE_MULT = "请确保Talbot倍数为正值。"
VAL_POSITIVE_ITER = "请确保计算次数(迭代次数)为正值。"

# ==================== Batch TLRec ====================
BATCH_TITLE = "批量重建"
BATCH_INFO = (
    "该工具可自动重建多次采集数据，无需手动加载到TL重建界面。\n\n"
    "选定的采集文件夹内需要以下结构：\n"
    "  - Acquisitions/AcquisitionX/Object/ -> 物体TIFF图像栈\n"
    "  - Acquisitions/AcquisitionX/Reference/ -> 可选的局部参考图像栈\n"
    "  - Acquisitions/Reference/ -> 可选的全局参考图像栈\n\n"
    "对于每个AcquisitionX，程序优先使用局部Reference/文件夹；\n"
    "否则回退到全局Acquisitions/Reference/。\n"
    "结果(DPC、相位、透射、暗场)保存在：\n"
    "  Acquisitions/AcquisitionX/Retrieved/"
)
BATCH_ACQ_FOLDER = "采集文件夹:"
BATCH_REC_ALGO = "重建算法:"
BTN_RUN_BATCH = "运行批处理"
BTN_BROWSE = "浏览..."
BATCH_MISSING_FOLDER = "缺少文件夹"
BATCH_MISSING_MSG = "请选择采集文件夹。"
BATCH_INVALID_FOLDER = "文件夹无效"
BATCH_INVALID_MSG = "选择的路径不是文件夹。"
BATCH_FINISHED = "批量重建完成。"
BATCH_FINISHED_TITLE = "批处理完成"

# ==================== Error Dialog ====================
ERR_TITLE = "错误"
ERR_REC_MSG = "重建过程中发生错误:\n{}"
ERR_GENERAL_MSG = "发生错误:\n{}"

# ==================== Help Window ====================
HELP_TITLE = "XPCIpy - 帮助 / 用户指南"
HELP_TAB_SIM = "仿真指南"
HELP_TAB_REC = "重建指南"
HELP_TAB_BATCH = "批量重建指南"
HELP_TAB_ABOUT = "关于"

# ==================== Misc ====================
MISSING_REF_TITLE = "缺少参考图像"
MISSING_REF_MSG = "请先加载或拖放参考图像栈，再加载物体图像栈。"
INVALID_FILES = "文件无效"
INVALID_FILES_MSG = "拖放的文件不符合预期扩展名。"
DROP_ERROR = "拖放错误"
DROP_ERROR_MSG = "无法处理拖放的文件:\n{}"
RECONSTRUCTING = "正在重建..."
LOADING = "加载中"

# ==================== Tooltips (Inline) ====================
TT_N_PIXELS = "模拟波前的像素数 (n x n)。"
TT_PIXEL_SIZE = "波前网格的像素尺寸，单位为微米。"
TT_BEAM_SHAPE = "光束几何形状：'Plane' = 平行光束，'Conical' = 发散锥形光束。"
TT_SPECTRUM = "选择能谱文件或 'Monoenergetic' 表示单一能量。"
TT_DSO = "源到物体的距离，单位为厘米。"
TT_FWHM_SOURCE = "源的半高全宽 (FWHM)，单位为微米。"
TT_ENERGY = "光束能量，单位为 keV。用于与波长相关的计算。"
TT_OBJECT = "要模拟的物体类型：'Sphere'（球体）或 'Cylinder'（圆柱体）。"
TT_IMAGE = "要模拟的图像类型：'Ideal'（理想探测器）或 'Realistic'（含探测器效应）。"
TT_DET_PX = "探测器像素尺寸，单位为微米。"
TT_DET_RES = "探测器分辨率，以半高全宽 (FWHM) 表示，单位为微米。"
TT_RUN_INLINE = "使用指定参数启动在线相衬成像仿真。"

# ==================== Tooltips (Check TL) ====================
TT_DESIGN_ENERGY = "设置的设计能量，单位为 keV。"
TT_PERIOD = "光栅周期，单位为微米。"
TT_DC = "光栅的占空比 (DC)，定义为栅条宽度与光栅周期的比值。"
TT_MATERIAL = "光栅栅条的材料（仅用于\"自定义\"光栅类型）。"
TT_BAR_HEIGHT = "光栅栅条高度，单位为微米（仅用于\"自定义\"光栅类型）。"
TT_GRATING_TYPE = "光栅类型：'Custom' 允许自定义参数，'Phase pi/2' 和 'Phase pi' 是标准相位光栅。"
TT_MULTIPLES = "Talbot距离的倍数（最大距离）。"
TT_ITERATIONS = "计算的距离数。"
TT_RUN_CHECK = "使用指定参数启动 Talbot-Lau 效应检查仿真。"

# ==================== Tooltips (TL Simulation) ====================
TT_DSG1 = "源到第一光栅 (G1) 的距离，单位为厘米。"
TT_DOG1 = "物体到第一光栅 (G1) 的距离，单位为厘米。"
TT_G1_PERIOD = "第一光栅 (G1) 的周期，单位为微米。"
TT_G1_PHASE = "第一光栅 (G1) 引入的相移。"
TT_MOVABLE = "选择在相位步进仿真中移动哪个光栅 (G1 或 G2)。"
TT_NUM_STEPS = "相位步进过程中的离散步数。"
TT_STEP_LENGTH = "每步的长度，单位为微米。"
TT_DET_RES_PX = "探测器分辨率，以像素尺寸表示，单位为微米。"
TT_RUN_TL = "使用指定参数启动 Talbot-Lau 相衬成像仿真。"

# ==================== Plot Axis Labels ====================
AXIS_PHASE_STEPPING = "相位步进"
AXIS_MULTIPLES_TD = "Talbot距离的倍数"
AXIS_INTENSITY = "强度"

# ==================== Overlay / Progress ====================
OVERLAY_LOADING = "加载中"
OVERLAY_RECONSTRUCTING = "正在重建..."

# ==================== ZIP Export ====================
ZIP_INLINE_TITLE = "XPCIpy 在线仿真"
ZIP_TL_TITLE = "Talbot-Lau 相衬仿真结果"

# ==================== Toggle Button ====================
TOGGLE_CREATE_ZIP = "创建仿真数据ZIP文件"

# ==================== Info Windows ====================
WIN_CITE_TITLE = "如何引用 XPCIpy"
WIN_CITE_HEADER = "如何引用 XPCIpy"
WIN_LICENSE_TITLE = "许可证"

# ==================== Help Content ====================
HELP_SIM_TITLE = "仿真模块"
HELP_INLINE_TITLE = "在线仿真"
HELP_INLINE_DESC = (
    "该模块模拟基于传播的（在线）相衬成像。\n"
    "该仿真遵循菲涅耳公式。\n"
    "主要参数：\n"
)
HELP_INLINE_BULLETS = [
    " - n：网格大小。\n",
    " - 网格像素尺寸 (µm)。\n",
    " - DSO / DOD 距离 (cm)。\n",
    " - 源 FWHM。\n",
    " - 光束能谱（单色/多色）。\n",
    " - 物体选择。\n",
    " - 探测器参数（像素尺寸、分辨率）。\n",
]
HELP_TL_SIM_TITLE = "Talbot-Lau 仿真"
HELP_TL_SIM_DESC = (
    "模拟带有 G1 和 G2 光栅的 Talbot-Lau 干涉仪。\n"
    "通过相位步进方法生成物体和参考图像栈。\n"
    "参数包括：\n"
)
HELP_TL_SIM_BULLETS = [
    " - 设计能量。\n",
    " - G1 周期和相位 (π, π/2)。\n",
    " - 相位步数。\n",
    " - 步长 (µm)。\n",
    " - 源参数（FWHM、能谱）。\n",
    " - 物体选择。\n",
    " - 探测器参数（像素尺寸、分辨率）。\n",
]
HELP_CHECK_TL_TITLE = "检查 Talbot 效应"
HELP_CHECK_TL_DESC = "可视化不同光栅配置下的 Talbot 毯。\n"

HELP_REC_TITLE = "Talbot-Lau 重建 (TLRec)"
HELP_REC_WORKFLOW = "相位步进重建工作流程"
HELP_REC_DESC = "重建过程使用通过相位步进获得的图像栈，在每个像素处拟合强度-相位曲线，提取物体参数。\n"
HELP_REC_BULLETS = [
    " - 加载图像栈（参考和物体）。确保调制曲线恰好覆盖一个周期。\n",
    " - 拟合算法：重建涉及将调制曲线拟合为正弦函数。常用选项包括 FFT（用于均匀 N）或最小二乘法。\n",
    " - 重建通道：\n",
    "   * **衰减（透射）**：标准吸收图像 ($I_0$)。\n",
    "   * **微分相位对比 (DPC)**：测量光束的角偏转，可用于重建累积相位 ($\\phi$)。\n",
    "   * **暗场成像 (DF)**：测量亚像素散射，与物体内部未分辨的微结构有关。\n",
]
HELP_REC_TOOLS_TITLE = "交互工具和诊断"
HELP_REC_TOOLS_DESC = "这些工具有助于验证原始数据和重建的质量：\n"
HELP_REC_TOOLS_BULLETS = [
    " - **调制曲线图**：点击图像显示该像素的强度-相位步曲线。用于诊断信号质量和噪声。\n",
    " - **窗宽/窗位调节**：调整重建图像的对比度和亮度，以查看特定细节。\n",
]

HELP_BATCH_TITLE = "批量重建工具"
HELP_BATCH_OVERVIEW = "概述"
HELP_BATCH_OVERVIEW_DESC = (
    "批量重建工具允许您自动处理多次采集数据，无需手动加载到 TLRec GUI。"
    "适用于大数据集、重复测量或扫描工作流程。\n"
)
HELP_BATCH_FOLDER_TITLE = "所需文件夹结构"
HELP_BATCH_FOLDER_DESC = "在所选的 **Acquisitions** 文件夹内，每个数据集必须遵循以下结构：\n"
HELP_BATCH_FOLDER_BULLETS = [
    " - **Acquisitions/AcquisitionX/Object/** -> 包含物体 TIFF 图像栈。\n",
    " - **Acquisitions/AcquisitionX/Reference/** -> 可选的局部参考图像栈。\n",
    " - **Acquisitions/Acquisitions/Reference/** -> 当局部参考不可用时使用的可选全局参考。\n",
]
HELP_BATCH_FOLDER_NOTE = (
    "如果采集目录内存在局部 *Acquisitions/Reference/* 文件夹，则使用该文件夹。\n"
    "否则，批处理工具将回退到全局参考。\n"
)
HELP_BATCH_WHAT_TITLE = "工具功能"
HELP_BATCH_WHAT_DESC = "对于每个有效的采集文件夹：\n"
HELP_BATCH_WHAT_BULLETS = [
    " - 加载 **Object/** 中的第一个 TIFF 图像栈。\n",
    " - 选择参考：**局部**（如果可用），否则使用**全局**。\n",
    " - 使用 TLRec 参数运行选定的相位提取算法。\n",
    " - 对微分相位应用维纳滤波器。\n",
    " - 将重建通道保存在：\n",
    "     *Acquisitions/AcquisitionX/Retrieved/*\n",
]
HELP_BATCH_OUTPUT_TITLE = "保存的输出"
HELP_BATCH_OUTPUT_BULLETS = [
    " - **DPC.tif** — 微分相位对比\n",
    " - **Phase.tif** — 维纳滤波后的积分相位\n",
    " - **Transmission.tif** — 吸收图像\n",
    " - **DarkField.tif** — 暗场/散射信号\n",
]
HELP_BATCH_HOW_TITLE = "如何使用批处理 GUI"
HELP_BATCH_HOW_BULLETS = [
    " - 选择 **Acquisitions** 文件夹。\n",
    " - 选择重建算法。\n",
    " - 点击**运行批处理**。\n",
]
HELP_BATCH_HOW_NOTE = (
    "日志窗口将报告：检测到的采集数据、使用的参考类型，"
    "以及每个已处理数据集的输出文件夹。\n"
)
HELP_BATCH_NOTES_TITLE = "注意事项与建议"
HELP_BATCH_NOTES_BULLETS = [
    " - 仅支持多页 TIFF 图像栈。\n",
    " - 如果存在多个 TIFF 文件，将按字母顺序选择第一个。\n",
    " - 原始数据不会被修改；所有结果存储在 *Retrieved/* 下。\n",
    " - 缺少 *Object/* 文件夹的采集数据将被跳过。\n",
]

HELP_ABOUT_TITLE = "关于 XPCIpy GUI"
HELP_ABOUT_DESC = (
    "\nXPCIpy 是一个用于X射线相衬成像仿真和重建的工具包。\n"
    "更多信息请访问: https://doi.org/10.1364/OE.573918\n"
)
HELP_ABOUT_MODULES = "模块"
HELP_ABOUT_MODULES_DESC = (
    " PCSim：在线 + Talbot-Lau 仿真\n"
    " TLRec：使用多种算法进行相位提取\n"
)
HELP_ABOUT_AUTHOR = "作者"
HELP_ABOUT_PROJECT = "项目"
HELP_ABOUT_PROJECT_DESC = " 本工作属于 PREDICO 项目，由西班牙科学与创新部资助 (PID2021-123390OB-C22)。\n"

# ==================== Warnings ====================
WARN_GRATING_SAMPLING = "光栅采样警告"
WARN_PHASE_STEPPING = "相位步进警告"

# ==================== Misc Status ====================
ST_FINISHED_SUFFIX = " 完成！"
ST_SIM_FINISHED = "仿真完成！"
ST_PREPARING_FIT = "正在准备 {} 的重建..."
ST_READY = "就绪"

# ==================== PBI Phase Retrieval (PBIRec) ====================
TAB_PBI_REC = "PBI相位重建"

# Drop zones
PBI_DZ_IMG_TITLE = "强度图像"
PBI_DZ_IMG_DESC = "拖放归一化强度TIFF图像（多距离请拖入多个文件或多页TIFF）"
PBI_DZ_IMG_BTN = "浏览强度图像"
PBI_DZ_FLAT_TITLE = "平场图像（可选）"
PBI_DZ_FLAT_DESC = "拖放平场/空场TIFF图像"
PBI_DZ_FLAT_BTN = "浏览平场图像"
PBI_DZ_DARK_TITLE = "暗场图像（可选）"
PBI_DZ_DARK_DESC = "拖放暗场TIFF图像"
PBI_DZ_DARK_BTN = "浏览暗场图像"

# Labels
PBI_LBL_DIST_TABLE = "传播距离 DOD（cm），每行一个值或逗号分隔"
PBI_DIST_PLACEHOLDER = "例如:\n5.0\n10.0\n20.0"
PBI_LBL_ENERGY = "能量 (keV)"
PBI_LBL_PIXEL = "探测器像素尺寸 (µm)"
PBI_LBL_DSO = "源到物体距离 DSO (cm)"
PBI_LBL_BEAM = "光束类型"
PBI_LBL_ALGO = "重建算法"
PBI_LBL_DELTA_BETA = "δ/β 比 (仅 Paganin)"
PBI_LBL_ALPHA = "正则化参数 α"
PBI_LBL_PAD = "边缘填充像素数"

# Buttons
PBI_BTN_SAVE_PHASE = "保存相位图"
PBI_BTN_SAVE_ABS = "保存吸收/厚度图"

# Tooltips
PBI_TT_ENERGY = "X射线光子能量，单位 keV。"
PBI_TT_PIXEL = "探测器像素的物理尺寸，单位 µm。"
PBI_TT_DSO = "源到物体距离（锥形束时需要），单位 cm。平行束请设为 0。"
PBI_TT_BEAM = "光束几何类型：Conical（锥形束/微焦点）或 Plane（平行束/同步辐射）。"
PBI_TT_ALGO = (
    "相位重建算法：\n"
    "• TIE 多距离 — 多距离传输强度方程（≥2张图）\n"
    "• CTF 多距离 — 多距离对比传递函数（弱物体近似）\n"
    "• Paganin 单距离 — 单距离同质材料假设"
)
PBI_TT_DELTA_BETA = "材料折射率实部与虚部之比 δ/β（仅 Paganin 方法使用）。软组织约 1000，塑料约 500-2000。"
PBI_TT_ALPHA = "Tikhonov 正则化参数，防止频域除零。值越大越平滑但会损失细节。"
PBI_TT_PAD = "在图像边缘填充的像素数，用于减少边缘伪影。0 = 不填充。"

# Status / messages
PBI_ST_RUNNING = "正在运行 PBI 相位重建 ({})..."
PBI_ST_FINISHED = "PBI 相位重建完成！"
PBI_ST_ZIP_SAVED = "PBI 重建结果 ZIP 已保存。"
PBI_IMAGES_LOADED = "张强度图像已加载"
PBI_ERR_LOAD_IMG = "无法加载图像文件。"
PBI_ERR_NO_IMAGES = "请先加载至少一张强度图像。"
PBI_ERR_PARAMS = "请确保能量和像素尺寸为正值。"
PBI_ERR_NEED_1_DIST = "Paganin 方法需要至少输入 1 个传播距离 (DOD)。"
PBI_ERR_DIST_MISMATCH = "图像数量 ({}) 与距离数量 ({}) 不匹配。多距离方法要求二者一一对应。"
PBI_NO_RESULT_TITLE = "无结果"
PBI_NO_RESULT_MSG = "尚未运行重建，无结果可保存。"

# Plot titles
PBI_PLOT_PHASE = "相位 φ(x,y)"
PBI_PLOT_ABS = "吸收 / 投影厚度"

# ==================== Help: PBI Rec ====================
HELP_TAB_PBI_REC = "PBI相位重建指南"
HELP_PBI_REC_TITLE = "PBI (传播成像) 相位重建"
HELP_PBI_REC_OVERVIEW = "概述"
HELP_PBI_REC_OVERVIEW_DESC = (
    "本模块用于从基于传播的 (Inline / PBI) X射线成像实验数据中重建相位信息。\n"
    "支持多距离与单距离两类方法，适用于锥形束和平行束几何。\n"
)
HELP_PBI_REC_ALGOS_TITLE = "支持的算法"
HELP_PBI_REC_ALGOS_BULLETS = [
    " - **TIE 多距离**：基于传输强度方程 (Transport of Intensity Equation)，需 ≥2 张不同传播距离的归一化强度图。适用于混合材料样品。\n",
    " - **CTF 多距离**：基于对比传递函数 (Contrast Transfer Function) 的弱物体近似，需 ≥2 张图。对弱相位/弱吸收物体效果较好。\n",
    " - **Paganin 单距离**：假设同质材料，仅需 1 张图 + δ/β 比值。简单快速，适用于单一材料样品。\n",
]
HELP_PBI_REC_WORKFLOW_TITLE = "工作流程"
HELP_PBI_REC_WORKFLOW_BULLETS = [
    " 1. 加载强度图像（拖放或浏览；多距离请加载多个文件或多页TIFF）。\n",
    " 2. 可选：加载平场 (flat) 和暗场 (dark) 图像用于归一化。\n",
    " 3. 在距离表中输入每张图对应的 DOD（物体到探测器距离，cm），每行一个值。\n",
    " 4. 填写能量 (keV)、像素尺寸 (µm)、DSO (cm)、光束类型等参数。\n",
    " 5. 选择算法并点击「运行」。\n",
    " 6. 查看/保存相位图和吸收图。\n",
]
HELP_PBI_REC_PARAMS_TITLE = "参数说明"
HELP_PBI_REC_PARAMS_BULLETS = [
    " - **能量 (keV)**：X射线光子能量。\n",
    " - **像素尺寸 (µm)**：探测器物理像素大小。\n",
    " - **DSO (cm)**：源到物体距离；平行束请设为 0。\n",
    " - **δ/β 比**：材料折射率实部与虚部之比（仅 Paganin 使用）。\n",
    " - **α (正则化)**：Tikhonov 正则化参数，防止噪声放大。\n",
    " - **填充像素**：边缘零填充像素数，减少边缘伪影。\n",
]

# ==================== PBI Preview ====================
PBI_PREVIEW_TITLE = "图像预览"
PBI_PREVIEW_INPUT = "输入图像"
PBI_PREVIEW_FLAT = "平场图像"
PBI_PREVIEW_DARK = "暗场图像"
PBI_PREVIEW_PHASE = "相位图"
PBI_PREVIEW_ABS = "吸收/厚度图"
PBI_PREVIEW_NONE = "（无图像）"
PBI_PREVIEW_IMG_N = "强度图 #{}"
PBI_PREVIEW_SHAPE = "尺寸: {}×{}"
PBI_PREVIEW_RANGE = "值域: [{:.4g}, {:.4g}]"
PBI_PREVIEW_PIXEL_INFO = "像素 ({}, {}): {:.6g}"
PBI_PREVIEW_WL_TITLE = "窗宽/窗位"
PBI_PREVIEW_WINDOW = "窗宽 (W)"
PBI_PREVIEW_LEVEL = "窗位 (L)"
PBI_PREVIEW_AUTO_WL = "自动"
PBI_PREVIEW_TAB_INPUT = "输入"
PBI_PREVIEW_TAB_RESULT = "结果"
