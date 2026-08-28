import logging

# Define handlers
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
console_handler.setFormatter(
    logging.Formatter("%(name)s | %(levelname)s | %(message)s")
)

file_handler = logging.FileHandler("DEBUG.log", mode="w")
file_handler.setLevel(logging.DEBUG)
file_handler.setFormatter(logging.Formatter("%(asctime)s | %(name)s | %(message)s"))

# Define logger
logger = logging.getLogger("SRC")
logger.setLevel(logging.INFO)
if not logger.handlers:
    logger.addHandler(console_handler)
    logger.addHandler(file_handler)

fea_logger = logging.getLogger("SRC.FEA")
iga_logger = logging.getLogger("SRC.IGA")
