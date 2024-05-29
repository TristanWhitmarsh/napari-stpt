def parse_channel_input(input_string):
    channels = []
    # Remove spaces and split input by commas
    input_string = input_string.replace(" ", "")
    parts = input_string.split(',')

    for part in parts:
        if '-' in part:
            # Handle ranges
            start, end = part.split('-')
            try:
                start = int(start)
                end = int(end)
                channels.extend(range(start, end + 1))
            except ValueError:
                print(f"Invalid range: {part}")
        else:
            # Single channel
            try:
                channels.append(int(part))
            except ValueError:
                print(f"Invalid channel: {part}")

    return channels